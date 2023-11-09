\name{deffS}
\alias{deffS}
\title{
Spencer design effect for an estimated total from a \emph{pps} sample
}
\description{
Compute the Spencer design effect for single-stage samples selected with probability proportional to a measure of size.
}
\usage{
deffS(p, w, y)
}
\arguments{
  \item{p}{vector of 1-draw selection probabilities, i.e., the probability that each unit would be selected in a sample of size 1.}
  \item{w}{vector of inverses of selection probabilities for a sample}
  \item{y}{vector of the sample values of an analysis variable}
}
\details{
The Spencer design effect is the ratio of the variance of the \emph{pwr-}estimator of the total of \eqn{y}, assuming that a single-stage sample is selected with replacement, to the variance of the total estimated in \emph{srswr}. Varying selection probabilities can be used.
}
\value{
  numeric design effect
}

\references{
Park, I., and Lee, H. (2004). Design Effects for the Weighted Mean and Total Estimators under Complex Survey Sampling. \emph{Survey Methodology}, 30, 183-193.

Spencer, B. D. (2000). An Approximate Design Effect for Unequal Weighting When Measurements May Correlate With Selection Probabilities. \emph{Survey Methodology}, 26, 137-138.

Valliant, R., Dever, J., Kreuter, F. (2018, chap. 14). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{deff}}, \code{\link{deffCR}}, \code{\link{deffH}}, \code{\link{deffK}}
}
\examples{
set.seed(-500398777)
    # generate population using HMT function
pop.dat <- as.data.frame(HMT())
mos <- pop.dat$x
pop.dat$prbs.1d <- mos / sum(mos)
    # select pps sample
require(sampling)
n <- 80
pk <- pop.dat$prbs.1d
sam <- UPrandomsystematic(pk)
sam <- sam==1
sam.dat <- pop.dat[sam, ]
dsgn.wts <- 1/pk[sam]
deffS(p=sam.dat$prbs.1d, w=dsgn.wts, y=sam.dat$y)
}

\keyword{methods}
\keyword{survey} 