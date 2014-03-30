\name{CVcalc2}
\alias{CVcalc2}

\title{
Coefficient of variation of an estimated total in a 2-stage sample
}
\description{
Compute the coefficient of variation of an estimated total in a two-stage design. Primary sampling units (PSUs) can be selected either with probability proportional to size (\emph{pps}) or with equal probability. Elements are selected via simple random sampling (\emph{srs}).
}
\usage{
CVcalc2(V=NULL, m=NULL , nbar=NULL, k=1, delta=NULL, Bsq=NULL, Wsq=NULL)
}

\arguments{
  \item{V}{
 unit relvariance of analysis variable in the population
 }
  \item{m}{
  number of sample PSUs
 }
  \item{nbar}{
  number of sample elements per PSU
 }
  \item{k}{
  ratio of \eqn{B^2 + W^2} to \eqn{V}. Default value is 1.
 }
  \item{delta}{
  measure of homogeneity equal to \eqn{B^2/(B^2 + W^2)}
 }
  \item{Bsq}{
  unit relvariance of PSU totals, equal to population variance of totals divided by \eqn{\bar{t}_{U}^2} if PSUs are selected by simple random sampling; or, equal to \eqn{S_{U1(pwr)}^2} divided by \eqn{t_{U}^2} if PSUs are selected by \emph{ppswr}
 }
  \item{Wsq}{
  within PSU relvariance, equal to average element population variance divided by \eqn{\bar{y}_{U}^2}
 }

}
\details{
\code{CVcalc2} computes the coefficient of variation of an estimated total for a two-stage sample.  PSUs can be selected either with varying probabilities
  and with replacement or with equal probabilities and with replacement. Elements within PSUs are selected by simple random sampling.
  The \eqn{CV} formula is appropriate for approximating the relvariance of the probability-with-replacement (\emph{pwr})-estimator of a total when
  the same number of elements is selected within each sample PSU.
}

\value{
 Value of the coefficient of variation of an estimated total
}

\references{
Cochran, W.G. (1977, pp.308-310). \emph{Sampling Techniques}. New York: John Wiley & Sons.

Saerndal,C.E., Swensson, B., and Wretman, J. (1992). \emph{Model Assisted Survey Sampling}. New York: Springer.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.2.1). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.

}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{CVcalc3}}
}
\examples{
CVcalc2(V=1, m=20 , nbar=5, k=1, delta=0.05)
CVcalc2(V=10, m=20 , nbar=5, k=1, delta=0.5)
CVcalc2(V=2.5, m=20 , nbar=5, k=2, Bsq=1, Wsq=4)
}

\keyword{methods}
\keyword{survey}
