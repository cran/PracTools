\name{nhis}
\alias{nhis}
\docType{data}
\title{
    National Health Interview Survey: Demographic variables
}
\description{
Demographic variables from a U.S. national household survey
}

\usage{data(nhis)}
\format{
  A data frame with 3,911 observations on the following 16 variables.
  \describe{
    \item{\code{ID}}{Identification variable}
    \item{\code{stratum}}{Sample design stratum}
    \item{\code{psu}}{Primary sampling unit, numbered within each stratum (1,2)}
    \item{\code{svywt}}{survey weight}
    \item{\code{sex}}{Gender (1 = male; 2 = female)}
    \item{\code{age}}{Age, continuous}
    \item{\code{age_r}}{Recoded age
               (3 = 18-24 years;
               4 = 25-44 years;
               5 = 45-64 years;
               6 = 65-69 years;
               7 = 70-74 years;
               8 = 75 years and older)
    }
    \item{\code{hisp}}{Hispanic ethnicity
               (1 = Hispanic;
               2 = Non-Hispanic)
     }
    \item{\code{marital}}{Marital status
               (1 = Separated;
               2 = Divorced;
               3 = Married;
               4 = Single/never married;
               5 = Widowed;
               9 = Unknown marital status)
     }
    \item{\code{parents}}{Parent(s) of sample person present in the family
               (1 = Mother, no father;
               2 = Father, no mother;
               3 = Mother and father;
               4 = Neither mother nor father)
     }
    \item{\code{parents_r}}{Parent(s) of sample person present in the family recode (1 = Yes; 2 = No)}
    \item{\code{educ}}{Education
               (1 = 8th grade or less;
               2 = 9-12th grade, no high school diploma;
               3 = High school graduate;
               4 = General education development (GED) degree recipient;
               5 = Some college, no degree;
               6 = Associate's degree, technical or vocational;
               7 = Associate's degree, academic program;
               8 = Bachelor's degree (BA, BS, AB, BBA);
               9 = Master's, professional, or doctoral degree)
     }
    \item{\code{educ_r}}{Education recode
               (1 = High school, general education development degree (GED), or less;
               2 = Some college;
               3 = Bachelor's or associate's degree;
               4 = Master's degree & higher)
     }
    \item{\code{race}}{Race (1 = White; 2 = Black; 3 = Other)}
    \item{\code{resp}}{Respondent (0 = nonrespondent; 1 = respondent)}
  }
}
\details{
The National Health Interview Survey (NHIS) is used to monitor health conditions in the U.S.
Data are collected through personal household interviews. Only demographic variables are included in this subset which was collected in 2003. The \code{nhis} data set contains observations for 3,911 persons. The file contains only persons 18 years and older.
}
\source{
National Health Interview Survey of 2003 conducted by the U.S. National Center for Health Statistics.
\url{http://www.cdc.gov/nchs/nhis.htm}
}

\seealso{
\code{\link{nhis.large}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
data(nhis)
str(nhis)
table(nhis$sex,nhis$age_r)
}
\keyword{datasets}
