\name{deffK}
\alias{deffK}
\title{
Kish design effect
}
\description{
Compute the Kish design effect due to having unequal weights.
}
\usage{
deffK(w)
}
\arguments{
  \item{w}{vector of inverses of selection probabilities for a sample}
}
\details{
The Kish design effect is the ratio of the variance of an estimated mean in stratified simple random sampling without replacement (\emph{stsrswor}) to the variance of the estimated mean in \emph{srswor}, assuming that all stratum unit variances are equal. In that case, proportional allocation with equal weighting is optimal. \code{deffK} equals \eqn{1 + relvar(w)} where \eqn{relvar} is relvariance of the vector of survey weights. This measure is not appropriate in samples where unequal weighting is more efficient than equal weighting.
}
\value{
  numeric design effect
}

\references{
Kish, L. (1965). \emph{Survey Sampling}. New York: John Wiley & Sons.

Kish, L. (1992). Weighting for unequal Pi.  \emph{Journal of Official Statistics}, 8, 183-200.

Valliant, R., Dever, J., Kreuter, F. (2018, chap. 14). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{deff}}, \code{\link{deffCR}}, \code{\link{deffH}}, \code{\link{deffS}}
}
\examples{
data(nhis)
w <- nhis$svywt
deffK(w)
}

\keyword{methods}
\keyword{survey} 