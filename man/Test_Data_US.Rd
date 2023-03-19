\name{Test_Data_US}
\alias{Test_Data_US}
\docType{data}
\title{
    Accounting data for some US cities with latitude and longitude of the city centroids
}
\description{
    A list of US cities with their latitude and longitude centroids and other data
}
\usage{data(Test_Data_US)}
\format{
  A data frame with 381 cities with the following variables:
  \describe{
     \item{\code{ID}}{Sequential ID field}
     \item{\code{State}}{State name}
     \item{\code{City}}{City name}
	 \item{\code{Count}}{Number of records in city}
     \item{\code{Amount}}{Total dollar amount of records}
     \item{\code{lat}}{latitude of the city center}
     \item{\code{long}}{longitude of the city center}
     \item{\code{Y}}{Artificial analysis variable}
     }	
  }
\details{
This population has 381 US cities with the latitude and longitude of the city center. It is used to illustrate the use of the \code{GeoDistPSU} and \code{GeoDistMOS} functions.
}
\seealso{
\code{\link{GeoDistPSU}},\code{\link{GeoDistMOS}}
}

\examples{
data(Test_Data_US)
str(Test_Data_US)
}
\keyword{datasets}
