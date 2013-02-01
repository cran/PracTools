\name{labor}
\alias{labor}
\docType{data}
\title{
    Labor force population
}
\description{
    A clustered population of persons extracted from the September 1976 Current Population Survey (CPS)
}
\usage{data(labor)}
\format{
  A data frame with 478 persons on the following variables:
  \describe{
     \item{\code{h}}{stratum}
     \item{\code{cluster}}{cluster (or segment) number}
	   \item{\code{person}}{person number}
     \item{\code{age}}{age of person}
     \item{\code{agecat}}{age category
         (1 = 19 years and under;
         2 = 20-24;
         3 = 25-34;
         4 = 35-64;
         5 = 65 years and over)
     }
     \item{\code{race}}{race	
         (1 = non-Black;
         2 = Black)
     }
     \item{\code{sex}}{Gender (1=Male; 2=Female)}
     \item{\code{HoursPerWk}}{Usual number of hours worked per week}
     \item{\code{WklyWage}}{Usual amount of weekly wages (in 1976 U.S. dollars)}	
     \item{\code{y}}{An artificial variable generated to follow a model with a common mean. Persons in the same cluster are correlated.  Persons in different clusters are uncorrelated under the model.
     }	
}
}
\details{
This population is a clustered population of 478 persons extracted from the
September 1976 Current Population Survey (CPS) in the United States.
The clusters are compact geographic areas used as one of the stages of sampling in the
CPS and are typically composed of about 4 nearby households.
The elements within clusters for this illustrative population are individual persons.
}

\source{
Current Population Survey, \url{http://www.census.gov/cps/}
}

\examples{
data(labor)
str(labor)
table(labor$h)
hist(labor$WklyWage)
}
\keyword{datasets}
