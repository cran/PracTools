\name{MDarea.popA}
\alias{MDarea.popA}
\docType{data}
\title{
    Maryland area population
}
\description{
    An artificial population of census tracts, block groups, and persons
}
\usage{data(MDarea.popA)}
\format{
  A data frame with 343,398 persons on the following variables:
  \describe{
     \item{\code{PSU}}{A grouping of block groups (\code{BLKGROUP}) which has about 4290 persons}
     \item{\code{SSU}}{A grouping of block groups which has about 1010 persons}
	 \item{\code{TRACT}}{A geographic area defined by the Census Bureau.  Tracts generally have between 1,500 and 8,000 people but have a much wider range in Anne Arundel county.}
     \item{\code{BLKGROUP}}{Block group. A geographic area defined by the Census Bureau.  Block groups generally have between 600 and 3,000 people.}
     \item{\code{Hispanic}}{Hispanic ethnicity (1=Hispanic; 2=Non-Hispanic)}
     \item{\code{Gender}}{ Gender (1=Male; 2=Female)}
     \item{\code{Age}}{23 level age category
        (1 = Under 5 years;
        2 = 5 to 9 years;
        3 = 10 to 14 years;
        4 = 15 to 17 years;
        5 = 18 and 19 years;
        6 = 20 years;
        7 = 21 years;
        8 = 22 to 24 years;
        9 = 25 to 29 years;
        10 = 30 to 34 years;
        11 = 35 to 39 years;
        12 = 40 to 44 years;
        13 = 45 to 49 years;
        14 = 50 to 54 years;
        15 = 55 to 59 years;
        16 = 60 and 61 years;
        17 = 62 to 64 years;
        18 = 65 and 66 years;
        19 = 67 to 69 years;
        20 = 70 to 74 years;
        21 = 75 to 79 years;
        22 = 80 to 84 years;
        23 = 85 years and over)
     }
     \item{\code{person}}{Counter for person within tract/block group/Hispanic/Gender/Age combination}
     \item{\code{y1}}{Artificial continuous variable}
     \item{\code{y2}}{Artificial continuous variable}
     \item{\code{y3}}{Artificial continuous variable}
     \item{\code{ins.cov}}{Medical coverage
            (0 = person does not have medical insurance coverage;
            1 = person has medical insurance coverage)
     }
     \item{\code{hosp.stay}}{Overnight hospital stay
            (0 = person did not have an overnight hospital stay in last 12 months;
            1 = person did have an overnight hospital stay in last 12 months)
     }
  }
}
\details{
A dataset of 343,398 persons based on the 2000 decennial U.S. Census for Anne Arundel County
in the US state of Maryland. Person records were generated based on counts from the 2000 census.
Individual values for each person were generated using models.
Groupings to form the variables \code{PSU} and \code{SSU} were done after sorting the census file by tract and block group within tract.

Note that \code{MDarea.popA} is different from the dataset, \code{MDarea.pop}, that is used in the book by Valliant, Dever, and Kreuter (2018).  \code{MDarea.pop} is larger with 403,997 persons. \code{MDarea.popA} was created by taking an equal probability, systematic subsample from \code{MDarea.pop}. \code{MDarea.popA} does have the same numbers of \code{TRACT}s, \code{PSU}s, and \code{SSU}s as \code{MDarea.pop}. The smaller data set was created to meet the CRAN size limit on installed packages. The full population, \code{MDarea.pop}, can be downloaded from \code{https://umd.app.box.com/v/PracTools2ndEdition}.
}

\source{
2000 U.S. decennial census, \code{http://www.census.gov/main/www/cen2000.html}
}
\references{
Valliant, R., Dever, J., Kreuter, F. (2018). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}
\examples{
data(MDarea.popA)
str(MDarea.popA)
table(MDarea.popA$TRACT)
table(MDarea.popA$TRACT, MDarea.popA$Hispanic)
}
\keyword{datasets}
