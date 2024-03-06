\name{nProp}
\alias{nProp}

\title{
Compute simple random sample size for estimating a proportion
}
\description{
Compute the simple random sample size for estimating a proportion based on different precision requirements.
}

\usage{
nProp(CV0 = NULL, V0 = NULL, pU = NULL, N = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CV0}{
target value of coefficient of variation of the estimated proportion
}
  \item{V0}{
target value of variance of the estimated proportion
}
  \item{pU}{
population proportion
}
  \item{N}{
number of units in finite population
}
}
\details{
The precision requirement of \eqn{p_s} can be set based on either a target coefficient of variation,
 \eqn{CV_0}, or a target variance, \eqn{V_0}. In either case, a value of \eqn{p_U} must be supplied.
}

\value{
  numeric sample size
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
}
\examples{
# srs sample size so that CV of estimated proportion is 0.05
# assuming the population is large and pU=0.01
# Both examples below are equivalent
nProp(V0=0.0005^2, N=Inf, pU=0.01) #or
nProp(CV0=0.05, N=Inf, pU=0.01)

# srswor sample size so that half-width of 2-sided 95% CI is 0.005
nProp(V0=(0.005/1.96)^2, N=Inf, pU=0.01)
}

\keyword{methods}
\keyword{survey} 