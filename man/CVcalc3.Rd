\name{CVcalc3}
\alias{CVcalc3}

\title{
Coefficient of variation of an estimated total in a 3-stage sample
}
\description{
Compute the coefficient of variation of an estimated total in a three-stage design. Primary sampling units (PSUs) can be selected
either with probability proportional to size (\emph{pps}) or with equal probability. Secondary units and elements within SSUs
are selected via simple random sampling (\emph{srs}).
}
\usage{
CVcalc3(V=NULL, m=NULL , nbar=NULL, qbar=NULL, k1=1, k2=1, delta1=NULL, delta2=NULL, 
            Bsq=NULL, Wsq=NULL, W2sq=NULL, W3sq=NULL)
}

\arguments{
 \item{V}{
 unit relvariance of analysis variable in the population
 }
 \item{m}{
  number of sample PSUs
 }
 \item{nbar}{
  number of sample secondary units per PSU
 }
 \item{qbar}{
  number of sample elements per SSU
 }
 \item{k1}{
  ratio of \eqn{B^2 + W^2} to \eqn{V}. Default value is 1.
 }
 \item{k2}{
  ratio of \eqn{W_{2}^{2} + W_{3}^{2}} to \eqn{V}. Default value is 1.
 }
 \item{delta1}{
  measure of homogeneity between PSUs equal to \eqn{B^2/(B^2 + W^2)} if PSUs are selected by simple random sampling; or, equal to \eqn{S_{U1(pwr)}^2} divided by \eqn{t_{U}^2} if PSUs are selected by \emph{ppswr}
 }
 \item{delta2}{
  measure of homogeneity between SSUs within PSUs, equal to \eqn{W_{2}^{2}/(W_{2}^{2} + W_{3}^{2})} 
 }
  \item{Bsq}{
  unit relvariance of PSU totals, equal to population variance of totals divided by \eqn{\bar{t}_{U}^2}
 }
  \item{Wsq}{
  within PSU relvariance of elements
 }
  \item{W2sq}{
  unit SSU relvariance
 }
  \item{W3sq}{
  unit element relvariance
 }

}
\details{
\code{CVcalc3} computes the coefficient of variation of an estimated total for a three-stage sample.  PSUs can be selected either with varying probabilities
  and with replacement or with equal probabilities and with replacement. SSUs and elements within SSUs are selected by simple random sampling.
  The \eqn{CV} formula is appropriate for approximating the relvariance of the probability-with-replacement (\emph{pwr})-estimator of a total when
  the same number of SSUs is selected in each PSU and the same number of elements is selected within each sample SSU.
}

\value{
 Value of the coefficient of variation of an estimated total
}

\references{
Cochran, W.G. (1977, pp.308-310). \emph{Sampling Techniques}. New York: John Wiley & Sons.

Saerndal,C.E., Swensson, B., and Wretman, J. (1992). \emph{Model Assisted Survey Sampling}. New York: Springer.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.2.4). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.

}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{CVcalc3}}
}
\examples{
CVcalc3(V=1, m=20 , nbar=5, qbar=10, delta1=0.02, delta2=0.10)
CVcalc3(V=1, m=20 , nbar=5, qbar=10, delta1=0.02, delta2=0.09)
CVcalc3(V=2, m=20 , nbar=5, qbar=10, k1=5, k2=10, Bsq=1, Wsq=9, W2sq=2 , W3sq=18 )
}

\keyword{methods}
\keyword{survey}
