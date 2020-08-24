binomial.funs <- function (gamma = 1, standardize = TRUE, intercept = TRUE,
                           lambda = NULL, nlambda = 50, lambda.min.ratio = 1e-04, cv = FALSE,
                           cv.rule = c("min", "1se")) {

  if (!require("glmnet", quietly = TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  if (!require("plyr", quietly = TRUE)) {
    stop("Package plyr not installed (required here)!")
  }
  check.num.01(gamma)
  check.bool(standardize)
  check.bool(intercept)
  if (is.null(lambda) && (length(nlambda) != 1 || round(nlambda) !=
                          nlambda || nlambda < 1 || nlambda > 100)) {
    stop("nlambda must be an integer between 1 and 100")
  }
  if (!is.null(lambda) && (!is.numeric(lambda) || min(lambda) <=
                           0 || (order(lambda) != length(lambda):1))) {
    stop("lambda must be a decreasing sequence of positive numbers")
  }
  check.pos.num(lambda.min.ratio)
  check.bool(cv)
  cv.rule = match.arg(cv.rule)
  if (cv) {
    train.fun = function(x, y, out = NULL) {
      return(cv.glmnet(x, y, alpha = gamma, nlambda = nlambda, family = "binomial",
                       lambda.min.ratio = lambda.min.ratio, lambda = lambda,
                       standardize = standardize, intercept = intercept))
    }
    predict.fun = function(out, newx) {
      return(predict(out$glmnet.fit, newx, type = "response", s = ifelse(cv.rule ==
                                                                           "min", out$lambda.min, out$lambda.1se)))
    }
    active.fun = function(out) {
      b = coef(out$glmnet.fit, s = ifelse(cv.rule == "min",
                                          out$lambda.min, out$lambda.1se))
      if (intercept)
        b = b[-1]
      return(list(which(b != 0)))
    }
  }
  else {
    train.fun = function(x, y, out = NULL) {
      return(glmnet(x, y, alpha = gamma, nlambda = nlambda, family = "binomial",
                    lambda.min.ratio = lambda.min.ratio, lambda = lambda,
                    standardize = standardize, intercept = intercept))
    }
    predict.fun = function(out, newx) {
      if (is.null(lambda)) {
        lambda = exp(seq(log(max(out$lambda)), log(min(out$lambda)),
                         length = nlambda))
      }
      return(predict(out, newx, type = "response", s = lambda))
    }
    active.fun = function(out) {
      if (is.null(lambda)) {
        lambda = exp(seq(log(max(out$lambda)), log(min(out$lambda)),
                         length = nlambda))
      }
      b = coef(out$glmnet.fit, s = lambda)
      if (intercept)
        b = b[-1, ]
      return(alply(b, 2, .fun = function(v) which(v !=
                                                    0)))
    }
  }
  return(list(train.fun = train.fun, predict.fun = predict.fun,
              active.fun = active.fun))
}


poisson.funs <- function (gamma = 1, standardize = TRUE, intercept = TRUE,
                          lambda = NULL, nlambda = 50, lambda.min.ratio = 1e-04, cv = FALSE,
                          cv.rule = c("min", "1se"))
{
  if (!require("glmnet", quietly = TRUE)) {
    stop("Package glmnet not installed (required here)!")
  }
  if (!require("plyr", quietly = TRUE)) {
    stop("Package plyr not installed (required here)!")
  }
  check.num.01(gamma)
  check.bool(standardize)
  check.bool(intercept)
  if (is.null(lambda) && (length(nlambda) != 1 || round(nlambda) !=
                          nlambda || nlambda < 1 || nlambda > 100)) {
    stop("nlambda must be an integer between 1 and 100")
  }
  if (!is.null(lambda) && (!is.numeric(lambda) || min(lambda) <=
                           0 || (order(lambda) != length(lambda):1))) {
    stop("lambda must be a decreasing sequence of positive numbers")
  }
  check.pos.num(lambda.min.ratio)
  check.bool(cv)
  cv.rule = match.arg(cv.rule)
  if (cv) {
    train.fun = function(x, y, out = NULL) {
      return(cv.glmnet(x, y, alpha = gamma, nlambda = nlambda,
                       lambda.min.ratio = lambda.min.ratio, lambda = lambda,
                       standardize = standardize, intercept = intercept))
    }
    predict.fun = function(out, newx) {
      return(predict(out$glmnet.fit, newx, s = ifelse(cv.rule ==
                                                        "min", out$lambda.min, out$lambda.1se)))
    }
    active.fun = function(out) {
      b = coef(out$glmnet.fit, s = ifelse(cv.rule == "min",
                                          out$lambda.min, out$lambda.1se))
      if (intercept)
        b = b[-1]
      return(list(which(b != 0)))
    }
  }
  else {
    train.fun = function(x, y, out = NULL) {
      return(glmnet(x, y, alpha = gamma, nlambda = nlambda, family = "poisson",
                    lambda.min.ratio = lambda.min.ratio, lambda = lambda,
                    standardize = standardize, intercept = intercept))
    }
    predict.fun = function(out, newx) {
      if (is.null(lambda)) {
        lambda = exp(seq(log(max(out$lambda)), log(min(out$lambda)),
                         length = nlambda))
      }
      return(predict(out, newx, s = lambda))
    }
    active.fun = function(out) {
      if (is.null(lambda)) {
        lambda = exp(seq(log(max(out$lambda)), log(min(out$lambda)),
                         length = nlambda))
      }
      b = coef(out$glmnet.fit, s = lambda)
      if (intercept)
        b = b[-1, ]
      return(alply(b, 2, .fun = function(v) which(v !=
                                                    0)))
    }
  }
  return(list(train.fun = train.fun, predict.fun = predict.fun,
              active.fun = active.fun))
}


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
##### fitting algorithm part

check.pos.num <- function (a)
{
  if (is.null(a) || length(a) != 1 || !is.numeric(a) || a <
      0)
    stop(paste(deparse(substitute(a)), "must be a positive number"))
}
check.bool <- function (b)
{
  if (is.null(b) || length(b) != 1 || !is.logical(b))
    stop(paste(deparse(substitute(b)), "must be a Boolean"))
}
check.num.01 <- function(a){
  if (is.null(a) || length(a) != 1 || !is.numeric(a) || a <
      0 || a > 1)
    stop(paste(deparse(substitute(a)), "must be a number between 0 and 1"))
}

check.args <- function(x = NULL, y = NULL, x0 = NULL, alpha = NULL, train.fun = NULL,
                       predict.fun = NULL, mad.train.fun = NULL, mad.predict.fun = NULL,
                       special.fun = NULL){
  if (is.null(x) || !is.numeric(x))
    stop("x must be a numeric matrix")
  if (is.null(y) || !is.numeric(y))
    stop("y must be a numeric vector")
  if (nrow(x) != length(y))
    stop("nrow(x) and length(y) must match")
  if (is.null(x0) || !is.numeric(x0))
    stop("x0 must be a numeric matrix")
  if (ncol(x) != ncol(x0))
    stop("ncol(x) and ncol(x0) must match")
  check.num.01(alpha)
  if (is.null(train.fun) || !is.function(train.fun))
    stop("train.fun must be a function")
  if (is.null(predict.fun) || !is.function(predict.fun))
    stop("predict.fun must be a function")
  if (!is.null(mad.train.fun) && !is.function(mad.train.fun))
    stop("mad.train.fun must be a function")
  if (!is.null(mad.predict.fun) && !is.function(mad.predict.fun))
    stop("mad.predict.fun must be a function")
  if ((!is.null(mad.train.fun) && is.null(mad.predict.fun)) ||
      (is.null(mad.train.fun) && !is.null(mad.predict.fun)))
    stop("mad.train.fun and mad.predict.fun must both be provided")
  if (!is.null(special.fun) && !is.function(special.fun))
    stop("special.fun must be a function")
}
