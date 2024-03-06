\name{gamEst}
\alias{gamEst}
\title{
Estimate variance model parameter \eqn{\gamma}
}
\description{
Regresses a \emph{y} on a set of covariates X where \eqn{Var_M(y)=\sigma^2x^\gamma} and then
regresses the squared residuals on \eqn{log(x)} to estimate \eqn{\gamma}.
}
\usage{
gamEst(X1, x1, y1, v1)
}
\arguments{
  \item{X1}{
matrix of predictors in the linear model for \code{y1}
}
  \item{x1}{
vector of \emph{x}'s for individual units in the assumed specification of \eqn{Var_M(y)}
}
  \item{y1}{
vector of dependent variables for individual units
}
  \item{v1}{
vector proportional to \eqn{Var_M(y)}
}
}
\details{
The function \code{gamEst} estimates the power \eqn{\gamma} in a model where the variance
 of the errors is proportional to \eqn{x^\gamma} for some covariate x.
 Values of \eqn{\gamma} are typically in [0,2].  The function is iteratively called by \code{\link{gammaFit}}, which is normally the function that an analyst should use.
}
\value{
The estimate of \eqn{\gamma}.
}
\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}
\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\seealso{
\code{\link{gammaFit}}
}
\examples{
data(hospital)
x <- hospital$x
y <- hospital$y

X <- cbind(sqrt(x), x)
gamEst(X1 = X, x1 = x, y1 = y, v1 = x)
}
\keyword{models}

