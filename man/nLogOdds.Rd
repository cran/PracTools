\name{nLogOdds}
\alias{nLogOdds}

\title{
Calculate simple random sample size for estimating a proportion}
\description{
Calculate the simple random sample size for estimating a proportion using the log-odds transformation.
}

\usage{
    nLogOdds(moe.sw, e, alpha=0.05, pU, N=Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{moe.sw}{
 switch for setting desired margin of error
    (1 = CI half-width on the proportion;
    2 = CI half-width on a proportion divided by \code{pU})
}
  \item{e}{
desired margin of error
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
The function accepts five parameters, which are the same ones as accepted by \code{\link{nPropMoe}}.
The desired margin of error can be specified as the CI half-width on the proportion (\code{moe.sw=1})
or as the CI half-width as a proportion of the population value \code{pU} (\code{moe.sw=2}).
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
\code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}, \code{\link{nCont}}

}
\examples{
nLogOdds(moe.sw=1, e=0.05, alpha=0.05, pU=0.2, N=Inf)
nLogOdds(moe.sw=2, e=0.05, alpha=0.05, pU=0.2, N=Inf)
}

\keyword{methods}
\keyword{survey} 