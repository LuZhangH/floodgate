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


  fg.out = fg.inference.binary(mu_X = mu_X, mu_Xk = mu_Xk, Y = Y[i2,], S = S, M_n = ifelse(useMC == TRUE, M_n, NULL),
                               alevel = alevel, verbose = verbose, test = test, one.sided = one.sided)

  fg.out$cpu.time = (proc.time() - begin.fg)[3]
  return(fg.out = fg.out)
}

#' @export
fg.inference.binary <- function(mu_X, mu_Xk, Y, S, M_n = NULL, alevel = 0.05, verbose = TRUE, test = "z", one.sided = TRUE){
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
