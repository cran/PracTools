\name{hospital}
\alias{hospital}
\docType{data}
\title{
Hospital Data
}
\description{
The \code{hospital} data file is a national sample of short-stay hospitals with fewer than 1000 beds.
}
\usage{data(hospital)}
\format{
  A data frame with 393 observations on the following 2 variables.
  \describe{
    \item{\code{y}}{Number of patients discharged by the hospital in January 1968}
    \item{\code{x}}{Number of inpatient beds in the hospital}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The \code{hospital} data are from the National Hospital Discharge Survey conducted by the U.S. National Center for Health Statistics.
The survey collects characteristics of inpatients discharged from non-Federal short-stay hospitals in the United States.
This population is from the January 1968 survey and contains observations on 393 hospitals.

}
\source{
    National Center for Health Statistics Hospital Discharge Survey of 1968.
}
\references{
    Herson, J. (1976). An Investigation of Relative Efficiency of Least-Squares Prediction to Conventional Probability Sampling Plans.  \emph{Journal of the American Statistical Association}, 71, 700-703.

}
\examples{
data(hospital)
str(hospital)
}
\keyword{datasets}
