###########################################################
##
##  Part1: Conditional model Y|X
##
###########################################################

#' Generate response variables
#
#' This function generates response variables Y under a given conditional model.
#'
#' @param X a n by p matrix, containing all the covariates
#' @param beta a p-dimensional vector, true regression coefficients
#' @param Ydist the type of conditional model
#' @return A matrix of size n-by-1 containing the response variables
#'
#' @family gen
#'
#' @details
#' There are four choices of the conditional model:
#' use "gaussian" to specify a linear model with gaussian noises,
#' use "laplace" to specify a linear model with laplacian noises,
#' use "binom" to specifiy a logistic model, with Y taking values in \{-1,1\}
#' use "poisson" to specify a log-poisson model
#'
#' @references TBD
#'
#' @import stats
#' @export
gen.Y <- function(X, beta, Ydist = "gaussian"){

  n = nrow(X)
  if(Ydist == "gaussian")
    Y = X %*% beta + rnorm(n)
  if(Ydist == "binom"){
    Y = rbinom(n,1, exp(X%*%beta) / (1+exp(X%*%beta) ) )
    Y = 2 * Y - 1
  }
  if(Ydist == "poisson")
    Y = log(rpois(n, exp(X %*% beta))+1)
  if(Ydist == "laplace"){
    sign_Y <- 2 * rbinom(n, 1, 0.5) - 1
    Y <- sign_Y * rexp(n, 2) + X %*% beta
  }

  Y = matrix(Y, nrow = n, ncol = 1)
  return(Y)
}


###########################################################
##
##  Part2 : Model-X
##
###########################################################

#### obtain HMM parameters from a randomly generated HMM model
#' @export
HMM.parm <- function(p, K_st, M){
### p: number of nodes
### K_st: number of hidden states
### M: number of emission states
### pInit: K_st vector
### Q: (p-1) * K_st * K_st array
### pEmit: p * M * K_st array

  pInit = rep(1/K_st, K_st)
  Q = array(stats::runif((p-1)*K_st*K_st), c(p-1,K_st,K_st))
  for(j in 1:(p-1)){
    Q[j,,] = Q[j,,] / rowSums(Q[j,,])
  }
  ### each Q[j,k1,] should be summed to one
  pEmit = array(stats::runif(p*M*K_st), c(p,M,K_st))
  for(j in 1:p) {
    pEmit[j,,] = t(t(pEmit[j,,]) / colSums(pEmit[j,,]))
  }
  ### pEmit[j,,k1] should be summed to one

  return(list(pInit = pInit, Q = Q, pEmit = pEmit))
}


#### obtain HMM parameters from a real genotype model
#' @importFrom SNPknock loadHMM
#' @export
Genotypes.parm <- function(p){
  ### pInit: K_st vector
  ### Q: (p-1) * K_st * K_st array
  ### pEmit: p * M * K_st array
  if(p > 2462) return(NULL)
  if(requireNamespace("SNPknock", quietly = TRUE)){
    r_file = system.file("extdata", "genotypes_rhat.txt", package = "SNPknock")
    alpha_file = system.file("extdata", "genotypes_alphahat.txt", package = "SNPknock")
    theta_file = system.file("extdata", "genotypes_thetahat.txt", package = "SNPknock")
    char_file  = system.file("extdata", "genotypes_origchars", package = "SNPknock")
    hmm = loadHMM(r_file, alpha_file, theta_file, char_file, compact=FALSE)

    pInit = hmm$pInit
    Q = hmm$Q[1:(p-1),,]
    pEmit = hmm$pEmit[1:p,,]
  } else {
    warning("Please install the SNPknock package.")
  }
  return(list(pInit = pInit, Q = Q, pEmit = pEmit))
}
