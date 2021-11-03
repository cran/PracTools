\name{nhispart}
\alias{nhispart}
\docType{data}
\title{
    National Health Interview Survey data from 2003: socioeconomic variables
}
\description{
Socioeconomic variables from a U.S. national household survey
}
\usage{data(nhispart)}
\format{
  A data frame with 3,924 observations on the following variables.
  \describe{
    \item{\code{HHX}}{Household identification variable}
    \item{\code{PX}}{Person identifier within household}
    \item{\code{STRATUM}}{Sample design stratum}
    \item{\code{PSU}}{Primary sampling unit, numbered within each stratum (1,2)}
    \item{\code{WTFA}}{survey weight}
    \item{\code{SEX}}{Gender (1 = male; 2 = female)}
    \item{\code{AGE_P}}{Age of persons; values are 18-85 (85 includes age 85 and older)}
    \item{\code{R_AGE1}}{Age group
                  (3 = 18-24 years;
                   4 = 25-44 years;
                   5 = 45-64 years;
                   6 = 65-69 years;
                   7 = 70-74 years;
                   8 = 75 years and over)
    }
    \item{\code{ORIGIN_I}}{Hispanic ethnicity
                  (1 = Hispanic;
                   2 = Non-Hispanic)
    }
    \item{\code{RACERPI2}}{Race grouped
                  (1 = White only;
                   2 = Black/African American only;
                   3 = American Indian or Alaska native (AIAN) only;
                   4 = Asian only;
                   5 = Race group not releasable;
                   6 = Multiple race)
    }
    \item{\code{MRACRPI2}}{Race detailed
                  (1  = White;
                   2  = Black/African American;
                   3  = Indian (American), Alaska Native;
                   9  = Asian Indian;
                   10 = Chinese;
                   11 = Filipino;
                   15 = Other Asian;
                   16 = Primary race not releasable;
                   17 = Multiple race, no primary race selected)
    }
    \item{\code{RACRECI2}}{White/Black
                   (1 = White;
                    2 Black;
                    3 All other race groups)
    }
    \item{\code{R_MARITL}}{Marital status
                   (1 = Married - spouse in household;
                    2 = Married - spouse not in household;
                    3 = Married - unknown whether spouse in household;
                    4 = Widowed;
                    5 = Divorced;
                    6 = Separated;
                    7 = Never married;
                    8 = Living with partner;
                    9 = Unknown marital status)
    }
    \item{\code{CDCMSTAT}}{CDC marital status
                    (1 = Mother, no father;
                     2 = Father, no mother;
                     3 = Mother and father;
                     4 = Neither mother nor father)
    }
    \item{\code{INCGRP}}{Total combined family income group
                    (1 = 0-$4999;
                     2 = $5000-$9999;
                     3 = $10000-$14999;
                     4 = $15000-$19999;
                     5 = $20000-$24999;
                     6 = $25000-$34999;
                     7 = $35000-$44999;
                     8 = $45000-$54999;
                     9 = $55000-$64999;
                     10 = $65000-$74999;
                     11 = $75000 and over;
                     12 = $20000 or more (no detail);
                     13 = Less than $20000 (no detail);
                     97 = Refused;
                     98 = Not ascertained;
                     99 = Don't know)
    }
    \item{\code{PARENTS}}{Parent(s) present in the family
                    (1 = Mother, no father;
                     2 = Father, no mother;
                     3 = Mother and father;
                     4 = Neither mother nor father)
    }
    \item{\code{EDUC_R1}}{Highest level of education attained
                    (1 = Less than high school graduate;
                     3 = High school graduate or general education development degree (GED);
                     5 = Some college, no degree;
                     6 = AA degree, technical or vocational or AA degree, academic program or Bachelor's degree (BA, BS, AB, BBA);
                     9 = Master's, professional, or doctoral degree)
    }
    \item{\code{RAT_CAT}}{Ratio of family income to poverty level
                    (1 = Under 0.50;
                     2 = 0.50 to 0.74;
                     3 = 0.75 to 0.99;
                     4 = 1.00 to 1.24;
                     5 = 1.25 to 1.49;
                     6 = 1.50 to 1.74;
                     7 = 1.75 to 1.99;
                     8 = 2.00 to 2.49;
                     9 = 2.50 to 2.99;
                     10 = 3.00 to 3.49;
                     11 = 3.50 to 3.99;
                     12 = 4.00 to 4.49;
                     13 = 4.50 to 4.99;
                     14 = 5.00 and over;
                     99 = Unknown)
    }
  }
}
\details{
The National Health Interview Survey (NHIS) is used to monitor health conditions in the U.S. Data are collected through personal household interviews. Socioeconomic variables are included in this subset along with household and person codes. The \code{nhispart} data set contains observations on 3,924 persons extracted from the 2003 U.S. survey. The file contains only persons 18 years and older.
}

\source{
National Health Interview Survey of 2003 conducted by the U.S. National Center for Health Statistics.
\url{https://www.cdc.gov/nchs/nhis.htm}
}


\examples{
data(nhispart)
str(nhispart)
table(nhispart$STRATUM, nhispart$PSU)
table(nhispart$RACERPI2, nhispart$RACRECI2, useNA="always")
}

\keyword{datasets}
