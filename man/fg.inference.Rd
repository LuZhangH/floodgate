% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/method.R
\name{fg.inference}
\alias{fg.inference}
\title{Core procedure of floodgate}
\usage{
fg.inference(
  S,
  mu_X,
  mu_Xk,
  Y,
  V_mean = NULL,
  useMC = TRUE,
  one.sided = TRUE,
  alevel = 0.05,
  test = "z",
  verbose = TRUE
)
}
\arguments{
\item{S}{a list of selected variables.}

\item{mu_X}{a list of kength |S|, whose element is the matrix of mu(X) with dimension n2-by-1.}

\item{mu_Xk}{a list of kength |S|, whose element is the matrix of mu(Xk) with dimension n2-by-K or n2-by-1.}

\item{Y}{a n2 by 1 matrix, containing the response variables of the inferene samples.}

\item{V_mean}{A vector of length |S|, whose element is the expected conditional variance term Var(Xj |X-j).}

\item{useMC}{whether to use Monte Carlo estimators of the conditional quantities (default: TRUE).=}

\item{one.sided}{whether to obtain LCB or p-values via the one-sided way (default: TRUE).}

\item{alevel}{confidence level (defaul: 0.05).}

\item{test}{type of the hypothesis test (defaul: "z").}

\item{verbose}{whether to show intermediate progress (default: FALSE).}
}
\value{
A list of three objects.
inf.out: a matrix of |S|-by-4, containing the p-values, LCI, UCI and the floodgate LCB for variable in S;
S: a list of selected variables.
}
\description{
This function produces floodgate LCBs for given fitted mu.
}
\references{
\insertRef{LZ-LJ:2020}{floodgate}
}
\seealso{
Other methods: 
\code{\link{calculate.V_mean}()},
\code{\link{calulate.mu_Xk}()},
\code{\link{fg.inference.binary}()},
\code{\link{fit.mu}()},
\code{\link{floodgate.binary}()},
\code{\link{floodgate}()},
\code{\link{inference_general}()}
}
\concept{methods}
