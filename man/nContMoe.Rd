\name{nContMoe}
\alias{nContMoe}
\title{
Compute a simple random sample size for an estimated mean of a continuous variable based on margin of error
}
\description{
Compute a simple random sample size using a margin of error specified as the half-width of a normal approximation confidence interval or the half-width relative to the population mean.
}
\usage{
nContMoe(moe.sw, e, alpha=0.05, CVpop=NULL, S2=NULL, ybarU=NULL, N=Inf)
}
\arguments{
  \item{moe.sw}{
switch for setting desired margin of error (1 = CI half-width on the mean;
2 = CI half-width on the mean divided by \eqn{\bar{y}_U})
}
  \item{e}{
desired margin of error; either \eqn{e=z_{1-\alpha/2}\sqrt{V(\bar{y}_s)}} or
\eqn{e=z_{1-\alpha/2}CV(\bar{y}_s)}
}
  \item{alpha}{1 - (confidence level)}
  \item{CVpop}{unit (population) coefficient of variation}
  \item{S2}{population variance of the target variable}
  \item{ybarU}{population mean of target variable}
  \item{N}{number of units in finite population}
}
\details{
If \code{moe.sw}=1, then \code{S2} must be provided.  If \code{moe.sw}=2, then either (i) \code{CVpop} or (ii) \code{S2} and \code{ybarU} must be provided.
}
\value{
  numeric sample size
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2013, chap. 3). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
}
\examples{
nContMoe(moe.sw=1, e=0.05, alpha=0.05, S2=2)
nContMoe(moe.sw=1, e=0.05, alpha=0.05, S2=2, N=200)
nContMoe(moe.sw=2, e=0.05, alpha=0.05, CVpop=2)
nContMoe(moe.sw=2, e=0.05, alpha=0.05, CVpop=2, N=200)
nContMoe(moe.sw=2, e=0.05, alpha=0.05, S2=4, ybarU=2)
}

\keyword{methods}
\keyword{survey} 