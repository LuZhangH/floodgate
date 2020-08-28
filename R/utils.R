#### calculate the mMSEgap value
#' Compute the mMSE gap value
#'
#' This function computes the mMSE gap value.
#'
#' @param beta a p-dimensional vector, true regression coefficients.
#' @param Xmodel model of the covaraites (default: "gaussian").
#' @param Ydist the type of conditional model.
#' @param sigma_X.list a list of length p, with each element being the variance of the conditional distribution.
#' @param X a n by p matrix, containing all the covariates.
#' @param nulls.list a list of length p, whose element is a (|i2|*K)-dimensional vector, which contains
#' K set of null samples.
#' @return A vector of p, with each element being the value of I for a given covariate.
#'
#' @family utils
#'
#' @export
compute.movi <- function(beta, Xmodel = "gaussian", Ydist = "gaussian",
                       sigma_X.list = NULL, X = NULL, nulls.list = NULL ){

  #### sqrt version for the target
  #### non-group version

  ## nulls.list: each element is (Kn)*Lj matrix
  p = length(beta)
  S_star = sort(which(beta!=0))
  J = length(S_star)
  nonpara.target = rep(0, p)

  if(Ydist %in% c("gaussian","laplace")){
    if(Xmodel == "gaussian" & !is.null(sigma_X.list))
      nonpara.target = abs(beta)*sqrt(unlist(sigma_X.list))
    else if(Xmodel %in% c("hmm","genotype", "haplotype") & !is.null(nulls.list) & !is.null(X)){

      n = nrow(X)
      K = (dim(nulls.list[[1]])[1])/n
      mu_Xk = vector(mode = "list", length = J)
      ##  each element of mu_Xk: n*K matrix
      for(j in 1:J){
        Gj = S_star[j]
        mu_Xk[[j]] = matrix(nulls.list[[Gj]],nrow = n, ncol = K, byrow = FALSE)
      }
      ## sigma_X similarly defined for hmm
      sigma_X_S_star = lapply(1:J, function(j)
        sum((mu_Xk[[j]] - rowMeans(mu_Xk[[j]]))^2)/n/(K-1)
      )

      nonpara.target[S_star] = abs(beta[S_star])*sqrt(unlist(sigma_X_S_star))
    }
  } else {

    n = nrow(X)
    K = (dim(nulls.list[[1]])[1])/n
    mu_Xk = vector(mode = "list", length = J)
    # if(Ydist == "poisson"){
    #   for(j in 1:J){
    #     Gj = S_star[j]
    #     common_part = as.vector(X[,-Gj]%*%beta[-Gj]) ## n-dim vector
    #     common_part = rep(common_part, K)  ## (Kn)-dim vector
    #     change_part = as.vector(nulls.list[[j]]*beta[Gj]) ## (Kn)-dim vector
    #     sum_part = common_part + change_part ## (Kn)-dim vector
    #     mu_Xk[[j]] = matrix(logpos.mean(exp(sum_part)),nrow = n, ncol = K)
    #   }
    # }
    if(Ydist == "binom"){
      for(j in 1:J){
        Gj = S_star[j]
        common_part = as.vector(X[,-Gj]%*%beta[-Gj]) ## n*1
        common_part = rep(common_part, K)
        change_part = as.vector(nulls.list[[j]]*beta[Gj])
        sum_part = common_part + change_part
        mu_Xk[[j]] = matrix(exp(sum_part)/(1+exp(sum_part)),nrow = n, ncol = K)
      }
    }
    ##  each element of mu_Xk: n*K matrix
    nonpara.target[S_star] = sapply(1:J, function(j)
      sum((mu_Xk[[j]] - rowMeans(mu_Xk[[j]]))^2)/n/(K-1))

    nonpara.target[S_star] = sqrt(nonpara.target[S_star])
  }

  return(nonpara.target = nonpara.target) ### p dimensional vector
}

# logpos.mean <- function(lambda){
#   ### lambda is a vector
#   t_seq = 0:(5*round(max(lambda)))
#   t_seq_expand = rep(t_seq,length(lambda))
#   lambda_expand = rep(lambda, each = length(t_seq))
#   tmp = log(1+t_seq_expand)*(lambda_expand^(t_seq_expand))/factorial(t_seq_expand)*exp(-lambda_expand)
#   val = colSums(matrix(tmp, ncol = length(lambda), byrow =  FALSE),na.rm = TRUE)
#   return(val)
# }


