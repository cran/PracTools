\name{nPropMoe}
\alias{nPropMoe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple random sample size for a proportion based on margin of error
}
\description{
Calculates a simple random sample size based on a specified margin of error.
}
\usage{
nPropMoe(moe.sw, e, alpha = 0.05, pU, N = Inf)
}
\arguments{
  \item{moe.sw}{
switch for setting desired margin of error (1 = CI half-width on the proportion;
2 = CI half-width on a proportion divided by \eqn{p_U})
}
  \item{e}{
%desired margin of error; either \eqn{e=z_{1-\alpha/2}\sqrt{V(\eqn{p_s})}} or
desired margin of error; either \eqn{e=z_{1-\alpha/2}\sqrt{V(p_s)}} or
\eqn{e=z_{1-\alpha/2}CV(p_s)}
}
  \item{alpha}{
1 - (confidence level)
}
  \item{pU}{
population proportion
}
  \item{N}{
number of units in finite population
}
}
\details{
The margin of error can be set as the half-width of a normal approximation
confidence interval, \eqn{e=z_{1-\alpha/2}\sqrt{V(p_s)}}, or as the half-width
of a normal approximation confidence interval divided by the population proportion,
\eqn{e=z_{1-\alpha/2}CV(p_s)}. The type of margin of error is selected by the
parameter \code{moe.sw} where \code{moe.sw=1} sets \eqn{e=z_{1-\alpha/2}\sqrt{V(p_s)}} and \code{moe.sw=2} sets i.e., \eqn{e=\frac{z_{1-\alpha/2}\sqrt{V(p_s)}}{p_U}}.
}
\value{
    numeric sample size
}
\references{
Valliant, R., Dever, J., Kreuter, F. (2013, chap. 3).  \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nWilson}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# srs sample size so that half-width of a 95% CI is 0.01
# population is large and population proportion is 0.04
nPropMoe(moe.sw=1, e=0.01, alpha=0.05, pU=0.04, N=Inf)

# srswor sample size for a range of margins of error defined as
# half-width of a 95% CI
nPropMoe(moe.sw=1, e=seq(0.01,0.08,0.01), alpha=0.05, pU=0.5)

# srswor sample size for a range of margins of error defined as
# the proportion that the half-width of a 95% CI is of pU
nPropMoe(moe.sw=2, e=seq(0.05,0.1,0.2), alpha=0.05, pU=0.5)
}
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{methods}
\keyword{survey} 