###########################################################
##
##  Sampling null variables
##
###########################################################

#### generate null variables for Gaussian X model
#' @import stats
#' @export
sample.gaussian.nulls <- function(X, S, K, gamma_X.list_S, sigma_X.list_S, verbose = FALSE){

  ## assume Xj|X-j ~ N( gamma_X %*% X-j,  sigma^2 )
  p = ncol(X)
  n = nrow(X)
  J = length(S)
  nulls.list_S = vector(mode = "list", length = J)
  ### if we are in the non-grouping case, nulls.list_S[[j]] contains a (Kn)*1 array
  ### if we are in the grouping case, nulls.list_S[[j]] contains a (Kn)*|Gj| array
  ### Gj is the j-th group
  ### (Kn) means: K=1,1:n; K=2,1:n; K=3,1:n .......
  if(length(unlist(S)) == length(S)){
    for(j in 1:J){
      Gj = S[[j]]
      if(verbose == TRUE){
        cat(sprintf("Generate Gaussian nulls for variable %i on %i samples with %i null replicates...\\n",
                    Gj, n, K),"\n")
      }
      Lj = length(Gj)
      gamma_X = gamma_X.list_S[[j]]  ## (p-1) * 1
      sigma_X = sigma_X.list_S[[j]]  ##  1-dim
      ### X[, -Gj] %*% gamma_X: n*(p-1) %*% (p-1) = n*1
      tmp = rnorm(n*K, 0 ,sqrt(sigma_X)) + as.vector(X[, -Gj] %*% gamma_X)
      ### (Kn)-dim vector + n-dim vector = (Kn)-dim vector
      nulls.list_S[[j]] = matrix(tmp, ncol=1)
      ### nulls.list_S[[j]]: (Kn)*1 matrix
      ### with patterns: K=1,1:n; K=2,1:n; K=3,1:n .......
    }
  } else{
    for(j in 1:J){
      Gj = S[[j]]
      Lj = length(Gj)
      gamma_X = gamma_X.list_S[[j]] ## (p - Lj)* Lj
      sigma_X = sigma_X.list_S[[j]] ##  Lj * Lj
      ### already vectorized, to be double checked
      proj_X = X[, -Gj] %*% gamma_X ## n*(p-Lj) %*% (p-Lj)*Lj = n * Lj
      ### (Kn)*Lj + (Kn)*Lj
      nulls.list_S[[j]] = MASS::mvrnorm(n*K, rep(0, Lj), sigma_X) + proj_X[rep(1:n, K),]
      ### nulls.list_S[[j]]: (Kn)*Lj matrix
    }
  }

  return(nulls.list_S = nulls.list_S)
}

#### generate null variables for Gaussian copula X model
#' @import stats
#' @export
sample.g_copula.nulls <- function(X, S, K, gamma_X.list_S, sigma_X.list_S, verbose = FALSE){
  ## for generating marginally uniform random variables from a multivariate gaussian distribution

  ## assume Xj|X-j ~ N( gamma_X %*% X-j,  sigma^2 )
  p = ncol(X)
  n = nrow(X)
  J = length(S)
  nulls.list_S = vector(mode = "list", length = J)
  ### if we are in the non-grouping case, nulls.list_S[[j]] contains a (Kn)*1 array
  ### if we are in the grouping case, nulls.list_S[[j]] contains a (Kn)*|Gj| array
  ### Gj is the j-th group
  ### (Kn) means: K=1,1:n; K=2,1:n; K=3,1:n .......
  if(length(unlist(S)) == length(S)){
    for(j in 1:J){
      Gj = S[[j]]
      if(verbose == TRUE){
        cat(sprintf("Generate Gaussian nulls for variable %i on %i samples with %i null replicates...\\n",
                    Gj, n, K),"\n")
      }
      Lj = length(Gj)
      gamma_X = gamma_X.list_S[[j]]  ## (p-1) * 1
      sigma_X = sigma_X.list_S[[j]]  ##  1-dim
      ### X[, -Gj] %*% gamma_X: n*(p-1) %*% (p-1) = n*1
      tmp = rnorm(n*K, 0 ,sqrt(sigma_X)) + as.vector(X[, -Gj] %*% gamma_X)
      ### (Kn)-dim vector + n-dim vector = (Kn)-dim vector
      tmp = 2*pnorm(tmp)-1
      nulls.list_S[[j]] = matrix(tmp, ncol=1)
      ### nulls.list_S[[j]]: (Kn)*1 matrix
      ### with patterns: K=1,1:n; K=2,1:n; K=3,1:n .......
    }
  } else {
    for(j in 1:J){
      Gj = S[[j]]
      Lj = length(Gj)
      gamma_X = gamma_X.list_S[[j]] ## (p - Lj)* Lj
      sigma_X = sigma_X.list_S[[j]] ##  Lj * Lj
      ### already vectorized, to be double checked
      proj_X = X[, -Gj] %*% gamma_X ## n*(p-Lj) %*% (p-Lj)*Lj = n * Lj
      ### (Kn)*Lj + (Kn)*Lj
      nulls.list_S[[j]] = MASS::mvrnorm(n*K, rep(0, Lj), sigma_X) + proj_X[rep(1:n, K),]
      ### nulls.list_S[[j]]: (Kn)*Lj matrix
    }
  }

  return(nulls.list_S = nulls.list_S)
}


#### generate null variables for Gaussian X model via co-sufficient sampling
#' @import stats
#' @export
cosuff.gaussian.nulls <- function(X, i2, n21, S, K, sigma_X.list_S, verbose = FALSE){

  ## assume Xj|X-j ~ N( gamma_X %*% X-j,  sigma^2 )
  p = ncol(X)
  n2 = length(i2)
  n22 = round(n2/n21)
  J = length(S)
  X = X[i2,]
  nulls.list_S = vector(mode = "list", length = J)
  cosuf.X_proj.list = vector(mode = "list", length = J)
  cosuf.sigma_X.list = vector(mode = "list", length = J)

  ### currently only implement the non-grouping case,
  ### nulls.list_S[[j]] contains a (Kn)*1 dimensional vector
  ### with (nulls.list_S[[j]])[1:n + (k-1)*n] being the k-th set of null samples
  if(length(unlist(S)) == length(S)){
    for(j in 1:J){
      Gj = S[[j]]
      ### Gj is the j-th group with |Gj|=1
      if(verbose == TRUE){
        cat(sprintf("Generate co-sufficient Gaussian nulls for variable %i on %i batches with batch size %i and %i null relicates \\n",
                    Gj, n21, n22, K),"\n")
      }
      Lj = length(Gj)
      nulls.list_S[[j]] = rep(0, n2*K)
      cosuf.X_proj.list[[j]] = matrix(0, nrow = n2, ncol = Lj)
      cosuf.sigma_X.list[[j]] = rep(0, n21)
      blocks = split(1:n2, ceiling(seq_along(1:n2)/n22))
      for(bat in 1:n21){
        U = cbind(rep(1,n22), X[blocks[[bat]], -Gj]) ## U: n22*p matrix
        H = U%*%solve(t(U)%*%U)%*%t(U) ## n22*n22 matrix
        proj_X = H%*%X[blocks[[bat]], Gj, drop = FALSE] ## (n22*n22)*(n22*1)
        sigma_X = as.vector(sigma_X.list_S[[j]])*(diag(n22) - H) ## (n22*n22) matrix
        entry.idx = rep(blocks[[bat]], K) + rep((1:K -1)*n2, each = n22) ## n22*K
        #cosuf.Sigma.chol = chol(sigma_X)
        ## directly we have: matrix(rnorm(n22*K),K,n22)%*% cosuf.Sigma.chol -> K*n22-> then transpose it
        tmp = as.vector((diag(n22) - H) %*% matrix(rnorm(n22*K, sd = sqrt(as.vector(sigma_X.list_S[[j]]))),n22, K)) + as.vector(proj_X) ## 2nd term: n22*1
        nulls.list_S[[j]][entry.idx] = matrix(tmp , ncol=1)
        cosuf.sigma_X.list[[j]][bat] = mean(diag(sigma_X)) ## diag(n22*n22)
        cosuf.X_proj.list[[j]][1:n22 + (bat-1)*n22,] = proj_X ## n22*1
      }
    }
    return(list(nulls.list_S = nulls.list_S, cosuf.X_proj.list = cosuf.X_proj.list, cosuf.sigma_X.list = cosuf.sigma_X.list))
  }
}


###################
#' @import stats
#' @export
cosuff.g_copula.nulls <- function(X, i2, n21, S, K, sigma_X.list_S, verbose = FALSE){

  ## assume Xj|X-j ~ N( gamma_X %*% X-j,  sigma^2 )
  p = ncol(X)
  n2 = length(i2)
  n22 = round(n2/n21)
  J = length(S)
  X = X[i2,]
  nulls.list_S = vector(mode = "list", length = J)
  ### currently only implement the non-grouping case,
  ### nulls.list_S[[j]] contains a (Kn)*1 dimensional vector
  ### with (nulls.list_S[[j]])[1:n + (k-1)*n] being the k-th set of null samples
  if(length(unlist(S)) == length(S)){
    for(j in 1:J){
      Gj = S[[j]]
      ### Gj is the j-th group with |Gj|=1
      if(verbose == TRUE){
        cat(sprintf("Generate co-sufficient Gaussian nulls for variable %i on %i batches with batch size %i and %i null relicates \\n",
                    Gj, n21, n22, K),"\n")
      }
      Lj = length(Gj)
      nulls.list_S[[j]] = rep(0, n2*K)
      # cosuf.X_proj.list[[j]] = matrix(0, nrow = n2, ncol = Lj)
      # cosuf.sigma_X.list[[j]] = rep(0, n21)
      blocks = split(1:n2, ceiling(seq_along(1:n2)/n22))
      for(bat in 1:n21){
        U = cbind(rep(1,n22), X[blocks[[bat]], -Gj]) ## U: n22*p matrix
        H = U%*%solve(t(U)%*%U)%*%t(U) ## n22*n22 matrix
        proj_X = H%*%X[blocks[[bat]], Gj, drop = FALSE] ## (n22*n22)*(n22*1)
        sigma_X = as.vector(sigma_X.list_S[[j]])*(diag(n22) - H) ## (n22*n22) matrix
        entry.idx = rep(blocks[[bat]], K) + rep((1:K -1)*n2, each = n22) ## n22*K
        #cosuf.Sigma.chol = chol(sigma_X)
        ## directly we have: matrix(rnorm(n22*K),K,n22)%*% cosuf.Sigma.chol -> K*n22-> then transpose it
        tmp = as.vector((diag(n22) - H) %*% matrix(rnorm(n22*K, sd = sqrt(as.vector(sigma_X.list_S[[j]]))),n22, K)) + as.vector(proj_X) ## 2nd term: n22*1
        tmp = 2*pnorm(tmp)-1
        nulls.list_S[[j]][entry.idx] = matrix(tmp , ncol=1)
        # cosuf.sigma_X.list[[j]][bat] = mean(diag(sigma_X)) ## diag(n22*n22)
        # cosuf.X_proj.list[[j]][1:n22 + (bat-1)*n22,] = proj_X ## n22*1
      }
    }

    return(list(nulls.list_S = nulls.list_S))
  }
}
