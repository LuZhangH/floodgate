#' Infer the model-free variable importance, MACM gap for binary responsess via floodgate
#
#' This function infers the model-free variable importance, MACM gap for binary responses via floodgate.
#'
#' @param X a n by p matrix, containing all the covariates.
#' @param Y a n by 1 matrix, containing the response variables.
#' @param i1 the index of training samples.
#' @param i2 the index of inference samples.
#' @param M_n the number of Monte Carlo samples for E(mu(X)|X-j).
#' @param nulls.list a list of length p, whose element is a (|i2|*(K))-dimensional vector, which contains
#' K set of null samples.
#' @param gamma_X.list a list of length p, with each element being the linear coefficient of
#' the given covariate on the other covariates (only relevant when Xmodel = "gaussian"; default: NULL).
#' @param sigma_X.list a list of length p, with each element being the variance of the conditional distribution.
#' @param Xmodel model of the covaraites (default: "gaussian").
#' @param funs a list of three: train.fun, active.fun and predict.fun.
#' @param algo a fitting algorithm (default: "lasso").
#' @param cv.rule indicates which rule should be used for the predict function, either "min" (the usual rule) or "1se" (the one-standard- error rule); default: "min").
#' See the glmnet help files for details.
#' @param one.sided whether to obtain LCB or p-values via the one-sided way (default: TRUE).
#' @param alevel confidence level (defaul: 0.05).
#' @param test type of the hypothesis test (defaul: "z").
#' @param verbose whether to show intermediate progress (default: FALSE).
#' @return A list of three objects.
#' inf.out: a matrix of |S|-by-4, containing the p-values, LCI, UCI and the floodgate LCB for variable in S;
#' S: a list of selected variables;
#' cpu.time: computing time.
#' @family methods
#'
#' @details TBD
#' @export
floodgate.binary <- function(X, Y, i1, i2, M_n = M_n,
                             nulls.list = NULL, gamma_X.list = NULL, sigma_X.list = NULL,
                             Xmodel = "gaussian", funs, algo = "lasso", cv.rule = "min",
                             one.sided = TRUE, alevel = 0.05, test = "z", verbose = TRUE){
  ### X: n*p Y: n vector
  ### i1: training sample index, i2: evaluation sample index

  ### Sigma: p * p matrix
  ### pInit: K_st vector
  ### Q: (p-1) * K_st * K_st array
  ### pEmit: p * M * K_st array

  begin.fg = proc.time()

  ### get the functions, given the algorithm
  train.fun = funs$train.fun
  active.fun = funs$active.fun
  predict.fun = funs$predict.fun

  ### fit mu_hat from the training samples
  n = nrow(X)
  p = ncol(X)
  n1 = length(i1)
  n2 = length(i2)

  model.fitted = fit.mu(X = X[i1,], Y = (Y[i1,]+1)/2, binary = TRUE, algo = algo,
                        train.fun = train.fun, active.fun = active.fun, verbose = verbose)

  out = model.fitted$out
  S = model.fitted$S
  S_unlist = unlist(S)

  if(length(S)==0){
    J = length(S)
    inf.out = matrix(0, J, 4)
    rownames(inf.out) = S
    colnames(inf.out) = c("P-value","LowConfPt", "UpConfPt","LCB")
    fg.out = list(inf.out = inf.out, S = S)
    fg.out$cpu.time = (proc.time() - begin.fg)[3]
    return(fg.out = fg.out)
  }

  ### calculate mu_Xk using nulls in general case
  ### much easider in special case: gaussian X and linear mu(x)
  if(Xmodel == "gaussian" & algo %in% c("lasso","ridge")){
    ### mu_Xk: n*1 matrix
    useMC = FALSE
    mu_Xk = calulate.mu_Xk(binary = TRUE, X = X, i2 = i2, S = S, out = out,
                           nulls.list_S = NULL, gamma_X.list_S = gamma_X.list[S_unlist],
                           useMC = useMC, Xmodel = Xmodel, algo = algo,
                           predict.fun = predict.fun, cv.rule = cv.rule, verbose = verbose)
  } else {
    ### mu_Xk: n*K matrix
    useMC = TRUE
    mu_Xk = calulate.mu_Xk(binary = TRUE, X = X, i2 = i2, S = S, out = out,
                           nulls.list_S = nulls.list[S_unlist], gamma_X.list_S = NULL,
                           useMC = useMC, Xmodel = Xmodel, algo = algo,
                           predict.fun = predict.fun, cv.rule = cv.rule, verbose = verbose)
  }

  ### fg core procedure
  if(algo %in% c("lasso","ridge")){
    mu_X = matrix(0, nrow = n2, ncol = 1)
  } else if(algo == "oracle"){
    mu_X = exp(X[i2,]%*%beta) / (1+exp(X[i2,]%*%beta))
    mu_X = 2*mu_X - 1
  } else {
    mu_X = matrix(predict.fun(out, X[i2,]), nrow = n2)
    mu_X = 2*mu_X - 1
  }
  ## mu_X: n2*1 matrix
  ## mu_Xk[[j]] : n2*K matrix or n2*1 matrix (when linear)


  fg.out = fg.inference.binary(S = S, mu_X = mu_X, mu_Xk = mu_Xk, Y = Y[i2,], M_n = ifelse(useMC == TRUE, M_n, NULL),
                               one.sided = one.sided, alevel = alevel, test = test, verbose = verbose)

  fg.out$cpu.time = (proc.time() - begin.fg)[3]
  return(fg.out = fg.out)
}

#' Core procedure of floodgate
#
#' This function produces floodgate LCBs for given fitted mu.
#'
#' @param S a list of selected variables.
#' @param mu_X a list of kength |S|, whose element is the matrix of mu(X) with dimension n2-by-1.
#' @param mu_Xk a list of kength |S|, whose element is the matrix of mu(Xk) with dimension n2-by-K or n2-by-1.
#' @param Y a n2 by 1 matrix, containing the response variables of the inferene samples.
#' @param M_n the number of Monte Carlo samples for E(mu(X)|X-j).
#' @param one.sided whether to obtain LCB or p-values via the one-sided way (default: TRUE).
#' @param alevel confidence level (defaul: 0.05).
#' @param test type of the hypothesis test (defaul: "z").
#' @param verbose whether to show intermediate progress (default: FALSE).
#' @return A list of three objects.
#' inf.out: a matrix of |S|-by-4, containing the p-values, LCI, UCI and the floodgate LCB for variable in S;
#' S: a list of selected variables.
#'
#' @family methods
#'
#' @details TBD
#' @export
fg.inference.binary <- function(S, mu_X, mu_Xk, Y, M_n = NULL, one.sided = TRUE, alevel = 0.05, test = "z", verbose = TRUE){
  ### Y: n2-dim vector
  ### mu_X: \mu evaluated at X,  n2*1
  ### mu_Xk: \mu evaluated at X_tilde then averaged over all K multiple samples
  ### J-list, mu_Xk[[j.idx]] contains a n2*K matrix
  ### S: set of selected variables/groups, S is a list, with j-th element Gj
  ### N: the normalization term

  ##
  if(!is.null(M_n) & !is.null(mu_Xk)){
    K = ncol(mu_Xk[[1]])
    if(K > M_n){
      tmp1.idx = 1:M_n
      tmp2.idx = setdiff(1:K, M_n)
    }
  }
  J = length(S)
  if(J > 0){
    inf.out = matrix(0, J, 4)
    rownames(inf.out) = S
    colnames(inf.out) = c("P-value","LowConfPt", "UpConfPt","LCB")
    if(verbose == TRUE){
      cat(sprintf("Performing fg analysis...\\n"),"\n")
    }
    if(is.null(M_n)){
      for(j in 1:J){
        R = 0.5 - (Y*(mu_X -  mu_Xk[[j]])<0)
        # mu_Xk[[j]]: n*K matrix, then rowMeans -> n*1
        inf.out[j,] = inference_general(R = 2*R, alevel = alevel, test = test, one.sided = one.sided)
      }

    } else{

      for(j in 1:J){
        mu_Xk_tmp1 = rowMeans(mu_Xk[[j]][,tmp1.idx])
        R = rowMeans((Y*(mu_Xk[[j]][,tmp2.idx] - mu_Xk_tmp1))<0) - (Y*(mu_X - mu_Xk_tmp1)<0)
        # mu_Xk[[j]]: n*K matrix, then rowMeans -> n*1
        inf.out[j,] = inference_general(R = 2*R, alevel = alevel, test = test, one.sided = one.sided)
      }

    }

    fg.out = list(inf.out = inf.out, S = S)
    class(fg.out) = "fg"
    return(fg.out = fg.out)
  }
}
