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
