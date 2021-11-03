\name{deffH}
\alias{deffH}
\title{
Henry design effect for \emph{pps} sampling and GREG estimation of totals
}
\description{
Compute the Henry design effect for single-stage samples when a general regression estimator is used for a total.
}
\usage{
deffH(w, y, x)
}
\arguments{
  \item{w}{vector of inverses of selection probabilities for a sample}
  \item{y}{vector of the sample values of an analysis variable}
  \item{x}{matrix of covariates used to construct a GREG estimator of the total of \eqn{y}.  This matrix does not include the intercept.}
}
\details{
The Henry design effect is the ratio of the variance of the general regression (GREG) estimator of a total of \eqn{y} to the variance of the estimated total in \emph{srswr}. Calculations for the Henry \emph{deff} are done as if the sample is selected in a single-stage and with replacement.  Varying selection probabilities can be used. The model for the GREG is assumed to be \eqn{y = \alpha + \beta x + \epsilon}, i.e., the model has an intercept.  }
\value{
  numeric design effect
}

\references{
Henry, K.A., and Valliant, R. (2015). A Design Effect Measure for Calibration Weighting in Single-stage Samples. \emph{Survey Methodology}, 41, 315-331.

Valliant, R., Dever, J., Kreuter, F. (2018, chap. 14). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{deff}}, \code{\link{deffCR}}, \code{\link{deffK}}, \code{\link{deffS}}
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
pk <- n * pop.dat$prbs.1d
sam <- UPrandomsystematic(pk)
sam <- sam==1
sam.dat <- pop.dat[sam, ]
dsgn.wts <- 1/pk[sam]
deffH(w=dsgn.wts, y=sam.dat$y, x=sam.dat$x)
}

\keyword{methods}
\keyword{survey} 