\name{nhis.large}
\alias{nhis.large}
\docType{data}
\title{
    National Health Interview Survey: Demographic and health variables
}
\description{
Demographic and health related variables from a U.S. national household survey
}
\usage{data(nhis.large)}
\format{
  A data frame with 21,588 observations on the following 18 variables.
  \describe{
    \item{\code{ID}}{Identification variable}
    \item{\code{stratum}}{Sample design stratum}
    \item{\code{psu}}{Primary sampling unit, numbered within each stratum (1,2)}
    \item{\code{svywt}}{survey weight}
    \item{\code{sex}}{Gender (1 = male; 2 = female)}
    \item{\code{age.grp}}{Age group
                  (1 = < 18 years;
                  2 = 18-24 years;
                  3 = 25-44 years;
                  4 = 45-64 years;
                  5 = 65+)
    }
    \item{\code{hisp}}{Hispanic ethnicity
                  (1 = Hispanic;
                  2 = Non-Hispanic White;
                  3 = Non-Hispanic Black;
                  4 = Non-Hispanic All other race groups)
    }
    \item{\code{parents}}{Parents present in the household
                  (1 = mother, father, or both present;
                  2 = neither present)
    }
    \item{\code{educ}}{Highest level of education attained
                  (1 = High school graduate, graduate equivalence degree, or less;
                  2 = Some college;
                  3 = Bachelor's or associate's degree;
                  4 = Master's degree or higher;
                  NA = missing)
}
    \item{\code{race}}{Race
                  (1 = White;
                  2 = Black;
                  3 = All other race groups)
}
    \item{\code{inc.grp}}{Family income group
                  (1 = < $20K;
                  2 = $20000-$24999;
                  3 = $25000-$34999;
                  4 = $35000-$44999;
                  5 = $45000-$54999;
                  6 = $55000-$64999;
                  7 = $65000-$74999;
                  8 = $75K+;
                  NA = missing)
}
    \item{\code{delay.med}}{Delayed medical care in last 12 months because of cost
                  (1 = Yes;
                  2 = No;
                  NA = missing)
}
    \item{\code{hosp.stay}}{Had an overnight hospital stay in last 12 months
                  (1 = Yes;
                  2 = No;
                  NA = missing)
}
    \item{\code{doc.visit}}{During 2 WEEKS before interview, did {person} see a doctor or
	other health care professional at a doctor's office, a clinic,
	an emergency room, or some other place? (excluding overnight hospital stay)?
                  (1 = Yes;
                  2 = No)
}
    \item{\code{medicaid}}{Covered by medicaid, a governmental subsidy program for the poor
                  (1 = Yes;
                  2 = No;
                  NA = missing)
}
    \item{\code{notcov}}{Not covered by any type of health insurance
                  (1 = Yes;
                  2 = No;
                  NA = missing)
}
    \item{\code{doing.lw}}{What was person doing last week?
                  (1 = Working for pay at a job or business;
                  2 = With a job or business but not at work;
                  3 = Looking for work;
                  4 = Working, but not for pay, at a job or business;
                  5 = Not working and not looking for work;
                  NA = missing)
}
    \item{\code{limited}}{Is the person limited in any way in any activities because of physical, mental or emotional problems?
                  (1 = Limited in some way;
                  2 = Not limited in any way;
                  NA = missing)
}
  }
}
\details{
The National Health Interview Survey (NHIS) is used to monitor health conditions in the U.S. Data are collected through personal household interviews. Demographic variables and a few health related variables are included in this subset. The \code{nhis.large} data set contains observations on 21,588 persons. \code{nhis.large} is a set of 21,588 persons extracted from the 2003 U.S. survey. The file contains only persons 18 years and older.
}

\source{
National Health Interview Survey of 2003 conducted by the U.S. National Center for Health Statistics.
\url{http://www.cdc.gov/nchs/nhis.htm}
}

\seealso{
\code{\link{nhis}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
data(nhis.large)
str(nhis.large)
table(nhis.large$stratum, nhis.large$psu)
table(nhis.large$delay.med, useNA="always")}

\keyword{datasets}
