\name{TPV}
\alias{TPV}
\docType{data}
\title{
TPV Data
}
\description{
\code{TPV} is an example data file for illustrating the use of certainty (take-all) units in sampling
}
\usage{data(TPV)}
\format{
  A data frame with 67 observations on the following 2 variables:
  \describe{
    \item{\code{Total.Pot.Value}}{a measure of size for each unit; for example, maximum potential amount spent on a contract, i.e. base price plus all options.}
    \item{\code{Y}}{an analytic variable for each unit}
  }
}
\details{
The \code{TPV} data are used as an example for \code{nContOpt} which determines the optimal split of a sample between take-all and non-take-all units.
}
\seealso{
\code{\link{nContOpt}}
}
\examples{
data(TPV)
str(TPV)
}
\keyword{datasets}
