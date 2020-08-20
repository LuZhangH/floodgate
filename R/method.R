#### fit a model to obtain the regression function estimator: mu
#' @export
fit.mu <- function(X, Y, binary = FALSE, algo = "lasso", train.fun, active.fun, verbose = TRUE){

  X = as.matrix(X)
  Y = as.matrix(Y)
  if(binary == TRUE & algo %in% c("rf","rf_binary")){
    Y = as.factor(Y)
  } else {
    Y = as.matrix(Y)
  }

  n = nrow(X)
  p = ncol(X)
  if(verbose == TRUE){
    cat(sprintf("Initial training on %i samples with %s algorithm...\\n",
                n, algo),"\n")
  }
  out = train.fun(X, Y)
  active = active.fun(out)
  S = sort(unique(unlist(active)), decreasing = FALSE)
  S = as.list(S)

  ### S will be a list, with each element containing one selected variable
  return(list(out = out, S = S))
}

#### calculate mu(Xk) using the null replicates
#### Xk: a set of null samples
#### mu: the regression function estimator
#' @export
calulate.mu_Xk <- function(binary = FALSE, X, i2, S, out, nulls.list_S = NULL, gamma_X.list_S = NULL,
                           useMC = TRUE, Xmodel = "gaussian", algo = "lasso",
                           predict.fun, cv.rule = "min", verbose = TRUE){

  ## nulls.list_S: a (Kn)*Lj matrix for each element in nulls.list_S
  J = length(S)
  n = nrow(X)
  n2 = length(i2)
  if(!is.null(nulls.list_S)){
    K = (dim(nulls.list_S[[1]])[1])/n
    infer_nulls.idx = rep(i2,K) + rep(n*(1:K - 1), each = n2)
  }

  mu_Xk = vector(mode = "list", length = J)
  if(useMC == FALSE){
    if(algo %in% c("lasso","ridge") & Xmodel == "gaussian" & !is.null(gamma_X.list_S)){
      beta_hat = coef(out, s = ifelse(cv.rule == "min", out$lambda.min,
                                      out$lambda.1se))[-1]
      for(j in 1:J){
        Gj = S[[j]]
        if(verbose == TRUE){
          cat(sprintf("Calculating mu(Xk) for variable %i with %s algorithm on %i samples without useMCimation...\\n",
                      Gj, algo, n2),"\n")
        }
        X_proj = X[i2,-Gj]%*% gamma_X.list_S[[j]]
        ### X[i2, -Gj] %*% gamma_X: n2*(p-Lj) %*% (p-Lj)*Lj = n2*Lj
        ### X_proj: n2*Lj matrix vector
        tmp = (X_proj - X[i2,Gj])%*%sign(beta_hat[Gj])
        ## since I assign mu(X) to be zero, hence it shoud be
        ## sign(beta_j)* ( E(Xj|X-j) - Xj)
        mu_Xk[[j]] = matrix(tmp, nrow = n2, ncol = 1, byrow = FALSE)
        ### n2*1 matrix

        ### X[i2, -Gj] %*% gamma_X: n2*(p-Lj) %*% (p-Lj)*Lj = n2*Lj
        ### change_part: n2*Lj %*% Lj*1 = n2-dim vector
        #tmp = change_part + common_mu_X[[j]] ## n2-dim vector
        #mu_Xk[[j]] = matrix(tmp, nrow = n2, ncol = 1, byrow = FALSE)
        ### n2*1 matrix
      }
    } else {
      stop("Wrong input...")
    }
  } else if (useMC == TRUE & !is.null(nulls.list_S)){
    if(algo %in% c("lasso","ridge")){
      beta_hat = coef(out, s = ifelse(cv.rule == "min", out$lambda.min,
                                      out$lambda.1se))[-1]
      for(j in 1:J){
        Gj = S[[j]]
        if(verbose == TRUE){
          cat(sprintf("Calculating mu(Xk) for variable %i with %s algorithm on %i samples and %i null samples...\\n",
                      Gj, algo, n2, K),"\n")
        }
        Xnulls = (nulls.list_S[[j]])[infer_nulls.idx] ## (K n2)*Lj
        common_part = X[i2,Gj] %*% matrix(sign(beta_hat[Gj]), ncol= 1) ## n2* Lj matrix(sign(beta_hat[Gj]),ncol= 1)
        change_part = as.vector(Xnulls %*% matrix(sign(beta_hat[Gj]), ncol= 1)) ## (K n2) - dim vector
        tmp =  change_part - rep(common_part, K)
        ## since I assign mu(X) to be zero, hence it shoud be
        ## sign(beta_j)* ( E(Xj|X-j) - Xj)
        mu_Xk[[j]] = matrix(tmp, nrow = n2, ncol = K, byrow = FALSE)
        ### n2*K matrix
        if(binary == TRUE) mu_Xk[[j]] = 2*mu_Xk[[j]] - 1
      }

    } else if(algo %in% c("sam", "samLL")){
      for(j in 1:J){
        mu_Xk[[j]] = matrix(0, nrow = n2, ncol = K, byrow = FALSE)
        Gj = S[[j]]
        if(verbose == TRUE){
          cat(sprintf("Calculating mu(Xk) for variable %i with %s algorithm on %i samples and %i null samples...\\n",
                      Gj,algo, n2, K),"\n")
        }
        # for(k in 1:K){
        #   Xk_expand = X[i2,] ## (Kn_2) * p matrix
        #   ## assign nulls.list_S[[j]] to the submatrix of Xk_expand
        #   ## nulls.list_S[[j]]: (Kn)*Lj matrix
        #   ## hence we also need to index with [rep(i2,K),]
        #   nulls_eval_idx_nonpar = i2 + rep(n*(k - 1), each = n2)
        #   Xk_expand[, Gj] = (nulls.list_S[[j]])[nulls_eval_idx_nonpar] ## (K n2)*Lj
        #   mu_Xk[[j]][,k] = matrix(predict.fun(out, Xk_expand), nrow = n2, ncol = 1, byrow = FALSE)
        #   ### n2*K matrix
        # }
        tmp = predict.fun(out = out, X = X[i2,],
                          Xkj = (nulls.list_S[[j]])[infer_nulls.idx], Gj = Gj)
        mu_Xk[[j]] = matrix(tmp, nrow = n2, ncol = K, byrow = FALSE)
        if(binary == TRUE) mu_Xk[[j]] = 2*mu_Xk[[j]] - 1
      }
    } else if(algo %in% c("binom_lasso","binom_ridge") & binary == TRUE){
        beta_hat = coef(out, s = ifelse(cv.rule == "min", out$lambda.min,
                                        out$lambda.1se))[-1]
        beta0_hat = coef(out, s = ifelse(cv.rule == "min", out$lambda.min,
                                         out$lambda.1se))[1]
        for(j in 1:J){
          Gj = S[[j]]
          if(verbose == TRUE){
            cat(sprintf("Calculating mu(Xk) for variable %i with %s algorithm on %i samples and %i null samples...\\n",
                        Gj,algo, n2, K),"\n")
          }
          common_part = as.vector(X[i2,-Gj]%*%beta_hat[-Gj]) ## n*1
          common_part = common_part + beta0_hat
          common_part = rep(common_part, K)
          change_part = as.vector(((nulls.list_S[[j]])[infer_nulls.idx])*beta_hat[Gj])
          sum_part = common_part + change_part
          mu_Xk[[j]] = matrix(2*exp(sum_part)/(1+exp(sum_part)) - 1, nrow = n2, ncol = K, byrow = FALSE)
          ### n2*K matrix
        }

    } else {
      for(j in 1:J){
        mu_Xk[[j]] = matrix(0, nrow = n2, ncol = K, byrow = FALSE)
        Gj = S[[j]]
        if(verbose == TRUE){
          cat(sprintf("Calculating mu(Xk) for variable %i with %s algorithm on %i samples and %i null samples...\\n",
                      Gj,algo, n2, K),"\n")
        }
        for(k in 1:K){
          Xk_expand = X[i2,] ## (Kn_2) * p matrix
          ## assign nulls.list_S[[j]] to the submatrix of Xk_expand
          ## nulls.list_S[[j]]: (Kn)*Lj matrix
          ## hence we also need to index with [rep(i2,K),]
          infer_nulls.subset.idx = i2 + rep(n*(k - 1), each = n2)
          Xk_expand[, Gj] = (nulls.list_S[[j]])[infer_nulls.subset.idx] ## (K n2)*Lj
          mu_Xk[[j]][,k] = matrix(predict.fun(out, Xk_expand), nrow = n2, ncol = 1, byrow = FALSE)
          ### n2*K matrix
        }
        if(binary == TRUE) mu_Xk[[j]] = 2*mu_Xk[[j]] - 1
     }

   }

  }

  return(mu_Xk = mu_Xk)  ## each element of mu_Xk: n2*K matrix or n2*1 matrix
}

### calculate the expected conditional variance term without using mu_Xk, denoted by V_mean
#' @export
calculate.V_mean <- function(S, algo = "lasso", cv.rule = "min", out = NULL,
                             Xmodel = "gaussian", sigma_X.list_S = NULL, verbose = TRUE){

  J = length(S)
  V_mean = rep(0, J)

  if(Xmodel == "gaussian" & algo %in% c("lasso", "ridge") & !is.null(out) & !is.null(sigma_X.list_S)){
    beta_hat = coef(out, s = ifelse(cv.rule == "min", out$lambda.min,
                                    out$lambda.1se))[-1]
    if(verbose == TRUE){
      cat(sprintf("Calculating V_mean for %s algorithm without Monte Carlo samples\\n",
                  algo),"\n")
    }

    if(length(unlist(S)) == length(S)){
      V_mean =  sign(beta_hat[unlist(S)])^2 * unlist(sigma_X.list_S)
    } else{
      V_mean = sapply(1:J, function(j){
        betaGj_hat = beta_hat[S[[j]]];
        t(betaGj_hat) %*% sigma_X.list_S[[j]] %*% betaGj_hat
      })
    }

  } else {
    V_mean = NULL
  }

  return(V_mean = V_mean)
}



#### floodgate main function: general use,
#### e.g. only need to known the conditional distribution of Xj | X-j
floodgate <- function(X, Y, i1, i2, nulls.list = NULL, gamma_X.list = NULL, sigma_X.list = NULL,
                      Xmodel = "gaussian", funs, algo = "lasso", cv.rule = "min",
                      one.sided = TRUE, alevel = 0.05, test = "z", verbose = TRUE){
  ### X: n*p matrix
  ### Y: n-dim vector
  ### i1: training sample index, i2: inference sample index

  ### Sigma: p * p matrix
  ### pInit: K_st vector
  ### Q: (p-1) * K_st * K_st array
  ### pEmit: p * M * K_st array

  begin.fg = proc.time()

  ### fetch the functions, given the algorithm
  train.fun = funs$train.fun
  active.fun = funs$active.fun
  predict.fun = funs$predict.fun

  ### fit mu from the training samples
  n = nrow(X)
  p = ncol(X)
  n1 = length(i1)
  n2 = length(i2)
  model.fitted = fit.mu(X = X[i1,], Y = Y[i1,], algo = algo,
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
    # gamma_X.list_S = gamma_X.list[S_unlist]
    # sigma_X.list_S = sigma_X.list[S_unlist]

    ### will be the exact conditional mean quantity
    mu_Xk = calulate.mu_Xk(X = X, i2 = i2, S = S, out = out,
                           nulls.list_S = NULL, gamma_X.list_S = gamma_X.list[S_unlist],
                           useMC = useMC, Xmodel = Xmodel, algo = algo,
                           predict.fun = predict.fun, cv.rule = cv.rule, verbose = verbose)
    ### calculate the expected conditional variance term without using mu_Xk, denoted by V_mean
    V_mean = calculate.V_mean(S = S, algo = algo, cv.rule = cv.rule, out = out,
                              Xmodel = Xmodel, sigma_X.list_S = sigma_X.list[S_unlist],
                              verbose = verbose)
  } else {
    ### mu_Xk: n*K matrix
    useMC = TRUE
    # nulls.list_S = nulls_list[S_unlist]
    mu_Xk = calulate.mu_Xk(X = X, i2 = i2, S = S, out = out,
                           nulls.list_S = NULL, gamma_X.list_S = gamma_X.list[S_unlist],
                           useMC = useMC, Xmodel = Xmodel, algo = algo,
                           predict.fun = predict.fun, cv.rule = cv.rule, verbose = verbose)
    ## set V_mean = NULL
    V_mean = calculate.V_mean(S = S, algo = algo, cv.rule = cv.rule, out = out,
                              Xmodel = Xmodel, sigma_X.list_S = sigma_X.list[S_unlist],
                              verbose = verbose)
  }

  ### fg core procedure
  if(algo %in% c("lasso","ridge")){
    mu_X = matrix(0, nrow = n2, ncol = 1)
  } else{
    mu_X = matrix(predict.fun(out, X[i2,]), nrow = n2)
  }


  ## mu_X : n2*1 matrix
  ## mu_Xk[[j]] : n2*K matrix,
  fg.out = fg.inference(S = S, mu_X = mu_X, mu_Xk = mu_Xk, Y = Y[i2,], V_mean = V_mean, useMC = useMC,
                        alevel = alevel, verbose = verbose, test = test, one.sided = one.sided)
  fg.out$cpu.time = (proc.time() - begin.fg)[3]

  return(fg.out = fg.out)
}


#### create a general version, e.g. works when Pj is known only
#' @export
fg.inference <- function(S, mu_X, mu_Xk, Y, V_mean = NULL, useMC = TRUE,
                         alevel = 0.05, verbose = TRUE, test = "z", one.sided = TRUE){
  ### Y: n2-dim vector
  ### mu_X: \mu evaluated at X,  n2*1
  ### mu_Xk: \mu evaluated at X_tilde then averaged over all K multiple samples
  ### J-list, mu_Xk[[j.idx]] contains a n2*K matrix
  ### S: set of selected variables/groups, S is a list, with j-th element Gj
  ### V_mean: the denominator term
  J = length(S)
  if(J > 0){
    inf.out = matrix(0, J, 4)
    rownames(inf.out) = S
    colnames(inf.out) = c("P-value","LowConfPt", "UpConfPt","LCB")

    if(useMC == FALSE){
      if(verbose == TRUE){
        cat(sprintf("Performing floodgate without Monte Carlo samples...\\n"),"\n")
      }
      for(j in 1:J){
        R = Y *(as.vector(mu_X) - rowMeans(mu_Xk[[j]]))
        # mu_Xk[[j]]: n2*1 matrix

        R = as.vector(R) ## n2-dim vector

        inf.out[j,] = inference_general(R = R, V = NULL, V_mean = V_mean,
                                        alevel = alevel, test = test, one.sided = one.sided)
      }
    } else {
      if(verbose == TRUE){
        cat(sprintf("Performing floodgate based on Monte Carlo samples...\\n"),"\n")
      }
      K = ncol(mu_Xk[[1]])
      for(j in 1:J){
        R = Y *(as.vector(mu_X) - rowMeans(mu_Xk[[j]]))
        # mu_Xk[[j]]: n2*K matrix, then rowMeans -> n2*1
        V = rowSums((mu_Xk[[j]] - rowMeans(mu_Xk[[j]]) )^2)/(K-1)

        R = as.vector(R) ## n2-dim vector
        V = as.vector(V) ## n2-dim vector

        inf.out[j,] = inference_general(R = R, V = V, V_mean = NULL,
                                        alevel = alevel, test = test, one.sided = one.sided)
      }

    }

    fg.out = list(inf.out = inf.out, S = S)
    class(fg.out) = "floodgate"
    return(fg.out = fg.out)
  }
}


inference_general <- function(R, V = NULL, V_mean = NULL, alevel  = 0.05, test = "z", one.sided = TRUE){
  ### V: 1/(K-1)*sum(Xi_tilde_k - mean(mu(Xi_tilde)))^2 , i \in [n2]
  ### R: Yi* (mu(Xi) - mean(mu(Xi_tilde)))  , i \in [n2]

  requireNamespace("stats", quietly = TRUE)

  n = length(R)

  if(is.null(V) & is.null(V_mean)){
    zz = R
    mm = mean(zz)
    ss = sd(zz)
  } else if(!is.null(V_mean) & is.null(V)){
    zz = R/sqrt(V_mean)
    mm = mean(zz)
    ss = sd(zz)
  } else if(!is.null(V) & is.null(V_mean)) {
    R_bar = mean(R)
    V_bar = mean(V)

    Sig11 = var(R)
    Sig12 = cov(R, V)
    Sig22 = var(V)

    s2 = ( Sig22*(R_bar/V_bar)^2/4 + Sig11 - Sig12*R_bar/V_bar )/V_bar

    mm = R_bar/sqrt(V_bar)
    ss = sqrt(s2)
  }

  ##
  if(one.sided == TRUE & test == "z"){
    ## since we know if we plug in -mu, we get theta(\mu) change the sign
    ## then mm with its sign changed and ss still the same
    ## since we already know theta(\mu_star) is nonnegative, hence we can just take
    ## the maximum of the inference results, done by using abs(mm)
    qq = qnorm(1 - alevel)
    lci = mm - qq * ss/sqrt(n)
    uci = mm + qq * ss/sqrt(n)
    #lcb =  abs(mm) - qq * ss/sqrt(n)
    lcb = max(0, lci) ## invert from one-sided tests
    pval = 1 - pnorm(mm/ss * sqrt(n))  ## one-sided pvalues
  } else if(one.sided == FALSE & test == "z"){
    ## since we know if we plug in -mu, we get theta(\mu) change the sign
    ## then mm with its sign changed and ss still the same
    ## since we already know theta(\mu_star) is nonnegative, hence we can just take
    ## the maximum of the inference results, done by using abs(mm)
    qq = qnorm(1 - alevel/2)
    lci = mm - qq * ss/sqrt(n)
    uci = mm + qq * ss/sqrt(n)
    #lcb =  abs(mm) - qq * ss/sqrt(n)
    lcb = max(lci,-uci,0) ## invert from two-sided tests
    pval = 1 - pnorm(abs(mm)/ss * sqrt(n)) ## two-sided pvalues
  }

  return(c(pval, lci, uci, lcb))
}
