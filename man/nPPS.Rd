\name{nPPS}
\alias{nPPS}
\title{
Calculate the sample size for a probability proportional to size (PPS) sample
}
\description{
Calculate the sample size for a probability proportional to size (PPS) sample, assuming the sample is selected with replacement.
}
\usage{
nPPS(X = NULL, Y = NULL, CV0 = NULL, V0 = NULL, N = NULL, V1 = NULL, ybarU = NULL)
}
\arguments{
  \item{X}{
    variable used for computing 1-draw probabilities; length is \eqn{N}, the population size; must be numeric
}
  \item{Y}{
    variable used for variance calculation; length is \eqn{N}, the population size; must be numeric
}
  \item{CV0}{
    target value of the coefficient of variation of the estimated total of \code{Y}
}
  \item{V0}{
    target value of the variance of the estimated total of \code{Y}; only one of \code{CV0} and \code{V0} can be specified
}
  \item{N}{
    population size; required if \code{X} or \code{Y} is NULL
}
  \item{V1}{
    unit variance for PPS calculation
}
  \item{ybarU}{
    population mean of \eqn{Y} (or an estimate of it)
}
}
\details{
    \code{nPPS} computes the sample size needed for a probability proportional to size sample or, more generally, a sample selected with varying probabilities, assuming the sample is selected with replacement (WR). Although these samples are rarely selected WR, the variance formula for WR samples is simple and convenient for sample size calculations. Population vectors can be input of \code{X}, a measure of size for selecting the sample, and \code{Y}, an analysis variable. Alternatively, the population size, \code{N}, the unit variance, \code{V1}, and the population mean of \eqn{Y}, \code{ybarU} can be inputs.
}
\value{
A list with four components:
    \item{N}{Size of the population}
    \item{V1}{Population variance of \eqn{Y} appropriate for a sample selected with varying probabilities and with replacement; see Valliant, Dever, and Kreuter (2018, sec. 3.4).}
    \item{ybarU}{Population mean of \eqn{Y}}
    \item{n}{Calculated sample size}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}
\author{
    George Zipf, Richard Valliant
}
\seealso{
\code{\link{nCont}}, \code{\link{nContMoe}}, \code{\link{nContOpt}}, \code{\link{unitVar}}
}
\examples{
library(PracTools)
data("smho.N874")
y <- smho.N874[,"EXPTOTAL"]
x <- smho.N874[, "BEDS"]
y <- y[x>0]
x <- x[x>0]
nPPS(X = x, Y = y, CV0 = 0.15)
nPPS(X = x, Y = y,  V0 = 2000000^2)
nPPS(CV0 = 0.15, N = length(y), V1 = (10^21), ybarU = mean(y))
}

\keyword{methods}
\keyword{survey} 