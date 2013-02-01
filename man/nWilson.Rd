\name{nWilson}
\alias{nWilson}

\title{
Calculate a simple random sample size for estimating a proportion
}
\description{
Calculate a simple random sample size for estimating a proportion using the Wilson method.
}
\usage{
nWilson(moe.sw, alpha = 0.05, pU, e)
}
\arguments{
  \item{moe.sw}{
 switch for setting desired margin of error (1 = CI half-width on the proportion;
 2 = CI half-width on a proportion divided by \code{pU})
}
  \item{alpha}{
 1 - (confidence level)
}
  \item{pU}{
 population proportion
}
  \item{e}{
 desired margin of error; either the value of CI half-width or the value of the half-width divided by \code{pU}
}
}
\details{
Calculate a simple random sample size using the Wilson (1927) method. A margin of error
can be set as the CI half-width on the proportion (\code{moe.sw=1}) or as the CI
half-width as a proportion of the population value \eqn{p_U} (\code{moe.sw=2}).
}
\value{
  \item{n.sam}{numeric sample size}
  \item{"CI lower limit"}{lower limit of Wilson confidence interval with computed sample size}
  \item{"CI upper limit"}{upper limit of Wilson confidence interval with computed sample size}
  \item{"length of CI"}{length of Wilson confidence interval with computed sample size}
}

\references{
Valliant, R., Dever, J., Kreuter,F. (2013, chap. 3).  \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.

Wilson, E.B. (1927). Probable inference, the law of succession, and statistical inference.  \emph{Journal of the American Statistical Association}, 22, 209-212.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}
}
\examples{
# srs sample size using Wilson method so that half-width of a 95% CI
# is 0.01. Population proportion is 0.04
nWilson(moe.sw = 1, pU = 0.04, e = 0.01)
}

\keyword{methods}
\keyword{survey} 