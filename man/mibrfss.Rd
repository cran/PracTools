\name{mibrfss}
\alias{mibrfss}
\docType{data}
\title{
    Michigan Behavioral Risk Factor Surveillance Survey
}
\description{
Demographic and health related variables from a U.S. household survey in the state of Michigan
}
\usage{data(mibrfss)}
\format{
  A data frame with 2485 observations on persons for the following 21 variables.
  \describe{
    \item{\code{SMOKE100}}{Smoked 100 or more cigarettes in lifetime
                  (1 = Yes;
                  2 = No)
}
    \item{\code{BMICAT3}}{Body mass index category
                  (1 = Neither overweight nor obese (BMI < 25);
                   2 = Overweight (25 <= BMI <= 30);
                   3 = Obese (BMI > 30) )
}
    \item{\code{AGECAT}}{Age group
                  (1 = 18-24 years;
                   2 = 25-34 years;
                   3 = 35-44 years;
                   4 = 45-54 years;
                   5 = 55-64 years;
                   6 = 65+)
}
    \item{\code{GENHLTH}}{General health (self-reported)
                  (1 = Excellent;
                   2 = Very good;
                   3 = Good;
                   4 = Fair;
                   5 = Poor)
}
    \item{\code{PHYSACT}}{Physical activity: In last month participated in activities such as running, calisthenics, golf, gardening, or walking for exercise
                  (1 = Yes;
                   2 = No)
}
    \item{\code{HIGHBP}}{High blood pressure: Have you ever been told by a doctor, nurse, or other health professional that you have high blood pressure?
                  (1 = Yes;
                   2 = No)
}
    \item{\code{ASTHMA}}{Asthma: Have you ever been told by a doctor, nurse, or other health professional that you have asthma?
                  (1 = Yes;
                   2 = No)
}
    \item{\code{HISPANIC}}{Hispanic ethnicity
                  (1 = Yes;
                  2 = No;
                  7 = Missing)
}
    \item{\code{WEIGHT}}{Body weight in pounds
}
    \item{\code{GENDER}}{Gender
                  (1 = Male;
                  2 = Female)
}
    \item{\code{CELLPHON}}{Has a wireless phone
                  (1 = Yes;
                  2 = No)
}
    \item{\code{INETHOME}}{Has access to the Internet at home
                  (1 = Yes;
                  2 = No)
}
    \item{\code{WEBUSE}}{How often do you use the Internet at home? Would you say, at least once a day, five to six times a week, two to four times a week, about once a week, less than once a week, or have you not used the Internet in the last month?
			     (1	= At least once a day;
                  2 = 5-6 times a week;
                  3	= 2-4 times a week;
                  4	= About once a week;
                  5	= Less than once a week;
                  6	= Not in the last month)
}
    \item{\code{RACECAT}}{Race
                  (1 = White;
                  2 = African American;
                  3 = Other)
}
    \item{\code{EDCAT}}{Education level
                  (1 = Did not graduate high school;
                  2 = Graduated high school;
                  3 = Attended college or technical school;
                  4 = Graduated from college or technical school)
}
    \item{\code{INCOMC3}}{Income category
                  (1 = Less than $15000;
                   2 = $15000 to less than $25000;
                   3 = $25000 to less than $35000;
                   4 = $35000 to less than $50000;
                   5 = $50000 or more)
}
    \item{\code{DIABETE2}}{Diabetes: Have you ever been told by a doctor, nurse, or other health professional that you have diabetes?
                  (1 = Yes;
                   2 = No)
}
    \item{\code{CHOLCHK}}{Cholesterol check: Blood cholesterol is a fatty substance found in the blood. Have you ever had your blood cholesterol checked?
                  (1 = Yes;
                   2 = No)
}
    \item{\code{BMI}}{Body mass index (continuous)
}
    \item{\code{BINGE2}}{Binge drinking: At risk for binge drinking based on alcohol consumption responses
                  (1 = Yes;
                   2 = No)
}
    \item{\code{ARTHRIT}}{Arthritis: Have you ever been told by a doctor, nurse, or other health professional that you have some form of arthritis, rheumatoid arthritis, gout, lupus, or fibromyalgia, or have joint symptoms of arthritis?
                  (1 = Yes;
                   2 = No;
                   3 = Don't know, not sure, or refused)
}
  }
}
\details{
The Michigan Behavioral Risk Factor Surveillance Survey (MIBRFSS) is part of a national state-by-state system of surveys used to monitor health conditions in the U.S. Data are collected through telephone household interviews. Demographic variables and a few health related variables are included in this subset. The \code{mibrfss} data set contains observations on 2845 persons and is  extracted from the 2003 U.S. survey. The file contains only persons 18 years and older.
}

\source{
Michigan Behavioral Risk Factor Surveillance Survey of 2003 sponsored by the U.S. Center for Disease Control.
\url{https://www.cdc.gov/brfss/}
}

\seealso{
\code{\link{nhis}}, \code{\link{nhis.large}}
}

\examples{
data(mibrfss)
str(mibrfss)
table(mibrfss$SMOKE100, useNA = "always")
table(mibrfss$BMICAT3, useNA="always")}

\keyword{datasets}
