this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)
#####################################################################################
#### specifying fitting algorithms
#### including lasso, step, ridge,
####           binom_lasso, binom_ridge, pois
####           rf, rf_binary,

#### loading the packages ####
library(methods)
library(conformalInference)
library(glmnet)
library(lars)
library(randomForest)
library(SAM)
source("algo_utils.R")

funs.list <- NULL

## Lasso, using cross-validation over 100 values of lambda, with the 1se rule
funs.list$lasso = lasso.funs(nlambda=100,cv=TRUE,cv.rule="min")

## Forward stepwise, using CV over 50 steps, with the 1se rule
funs.list$step = lars.funs(type="step",max.steps=50,cv=TRUE,cv.rule="min")

## ridge estimator
funs.list$ridge = elastic.funs(gamma=0,nlambda=100,cv=TRUE,cv.rule ="min" )

## logistic regression with lasso penalty
funs.list$binom_lasso = binomial.funs(nlambda=100,cv=TRUE,cv.rule="min")
funs.list$binom_ridge = binomial.funs(gamma = 0, nlambda=100,cv=TRUE,cv.rule="min")

## poisson regression with lasso penalty
funs.list$pois = poisson.funs(nlambda=100,cv=TRUE,cv.rule="min")


## random forests
my.rf.funs = rf.funs()
my.rf.funs$active.fun = function(out){
  return(list(sort(order(out$importance,decreasing=TRUE)[1:50])))
}
# We hijacked the active set function, take top 20 variables according to
# the forest important measure
funs.list$rf = my.rf.funs

my.rf.funs$predict.fun = function(out, newx){
  prob.mat = predict(out, newx, type = "prob")
  ## binary case: prob.mat has 2 columns: prob(Y=0|X) , prob(Y=1|X)
  return(prob.mat[,2])
  ## prob(Y=1|X)
}
funs.list$rf_binary = my.rf.funs

## Sparse additive model
my.predict.fix = function(object, newdata) {
  gcinfo(FALSE)
  out = list()
  nt = nrow(newdata)
  d = ncol(newdata)
  X.min.rep = matrix(rep(object$X.min, nt), nrow = nt, byrow = T)
  X.ran.rep = matrix(rep(object$X.ran, nt), nrow = nt, byrow = T)
  newdata = (newdata - X.min.rep)/X.ran.rep
  newdata = pmax(newdata, 0)
  newdata = pmin(newdata, 1)
  m = object$p * d
  Zt = matrix(0, nt, m)
  for (j in 1:d) {
    tmp = (j - 1) * object$p + c(1:object$p)
    unique.vals = unique(newdata[,j])
    tmp.mat = ns(unique.vals, df = object$p,
                 Boundary.knots = object$Boundary.knots[, j])
    for (i in 1:nrow(newdata)) {
      Zt[i, tmp] = tmp.mat[which(unique.vals==newdata[i,j]),]
    }
  }
  out$values = cbind(Zt, rep(1, nt)) %*% rbind(object$w, object$intercept)
  rm(Zt, newdata)
  return(out)
}

my.sam.train.fun = function(x, y, out=NULL, nfolds=10) {
  out = samQL(x,y)
  n = nrow(x)
  m = length(out$lambda)

  # Run CV manually
  ii = sample(rep(1:nfolds,length=n))
  cv.mat = matrix(0,n,m)
  for (k in 1:nfolds) {
    out.k = samQL(x[ii!=k,], y[ii!=k])
    pred.mat = my.predict.fix(out.k,x[ii==k,])$values
    cv.mat[ii==k,] = (pred.mat - y[ii==k])^2
  }
  cv.err = colMeans(cv.mat)
  cv.se = apply(cv.mat,2,sd)/sqrt(nfolds)

  # Compute i.min and i.1se, the min and 1se choices from CV
  i.min = which.min(cv.err)
  i.1se = min(which(cv.err <= cv.err[i.min] + cv.se[i.min]))

  # Append i.min and i.1se to the original SAM object
  out$i.min = i.min
  out$i.1se = i.1se
  return(out)
}

my.sam.active.fun = function(out, cv.rule="min") {
  i = ifelse(cv.rule=="min", out$i.min, out$i.1se)
  return(list(which(out$func[,i] != 0)))
}

fast.predict.samQL <- function(out, X, Xkj = NULL, Gj = 1, cv.rule = "min"){
  i = ifelse(cv.rule=="min", out$i.min, out$i.1se)
  coef_vec = rbind(out$w, out$intercept)
  coef_vec = coef_vec[,i, drop = FALSE]
  gcinfo(FALSE)
  d = ncol(X)
  nt = nrow(X)
  X.min.rep = matrix(rep(out$X.min, nt), nrow = nt, byrow = T)
  X.ran.rep = matrix(rep(out$X.ran, nt), nrow = nt, byrow = T)
  X = (X - X.min.rep)/X.ran.rep
  X = pmax(X, 0)
  X = pmin(X, 1)
  m = out$p * d
  Zt = matrix(0, nt, m)
  for (j in setdiff(1:d, Gj)) {
    tmp = (j - 1) * out$p + c(1:out$p)
    Zt[, tmp] = ns(X[, j], df = out$p, knots = out$knots[, j], Boundary.knots = out$Boundary.knots[, j])
  }
  out.mat = cbind(Zt, rep(1, nt)) %*% coef_vec
  j = Gj
  tmp = (j - 1) * out$p + c(1:out$p)

  if(is.null(Xkj)){
    Zt_tmp = ns(X[, j], df = out$p, knots = out$knots[, j], Boundary.knots = out$Boundary.knots[, j])
    out.mat = out.mat + Zt_tmp%*%coef_vec[tmp]
  } else {
    Xkj = (Xkj - X.min.rep[,j])/(X.ran.rep[,j])
    Xkj = pmax(Xkj,0)
    Xkj = pmin(Xkj,1)
    Zt_tmp = ns(Xkj, df = out$p, knots = out$knots[, j], Boundary.knots = out$Boundary.knots[, j])
    out.mat = rep(out.mat, K)
    out.mat = out.mat + Zt_tmp%*%coef_vec[tmp]
  }
  rm(Zt, X, Zt_tmp)
  return(out.mat)
}

fast.sam.funs = list(train.fun=my.sam.train.fun,
                     predict.fun=fast.predict.samQL,
                     active.fun=my.sam.active.fun)
funs.list$sam = fast.sam.funs

#### save them to a RData file

save(funs.list,
     file = "algos.RData")
