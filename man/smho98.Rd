\name{smho98}
\alias{smho98}
\docType{data}
\title{
Survey of Mental Health Organizations Data
}
\description{
Data from the 1998 Survey of Mental Health Organizations (SMHO)
}
\usage{data(smho98)}
\format{
  A data frame with 875 observations on the following variables.
  \describe{
    \item{\code{STRATUM}}{Sample design stratum
       (1 = Psychiatric Hospital, private;
       2 = Psychiatric Hospital, public;
       3 = Residential, children;
       4 = Residential, adults;
       5 = General Hospital, public, inpatient or residential care;
       6 = General Hospital, public, outpatient care only;
       7 = General Hospital, private, inpatient or residential care;
       8 = General Hospital, private, outpatient care only;       
       9 = Military Veterans, inpatient or residential care;
       10 = Military Veterans, outpatient care only;
       11 = Partial Care
       12 = Outpatient care, private;
       13 = Outpatient care, public;       
       14 = Multi-service, private;
       15 = Multi-service, public;       
       16 = Substance Abuse)
       }
    \item{\code{BEDS}}{Total inpatient beds}
    \item{\code{EXPTOTAL}}{Total expenditures in 1998}
    \item{\code{SEENCNT}}{Unduplicated client/patient count seen during year}
    \item{\code{EOYCNT}}{End of year count of patients on the role}
    \item{\code{Y_IP}}{Number of inpatient visits during year}
    \item{\code{OPCSFRST}}{Number of outpatients on the rolls on the first day of the reporting year}
    \item{\code{OPCSADDS}}{Number of outpatients admitted, readmitted, or transferred to the organization during the reporting year for less than a 24 hour period and not overnight}
    \item{\code{OPCSVIST}}{Number of outpatient visits during the reporting year for less than a 24 hour period and not overnight}
    \item{\code{EMGWALK}}{Number of emergency walk-ins during the reporting year}
    \item{\code{PSYREHAB}}{Number of visits for psychiatric rehabilitation services}
    \item{\code{IPCSADDS}}{Number of residential patients added during the reporting year or patients admitted for more than a 24 hour period}
  }
}
\details{
The 1998 SMHO was conducted by the U.S. Substance Abuse and Mental Health Services Administration.
It collected data on mental health care organizations and general hospitals that provide mental health care services,
with an objective to develop national and state level estimates for total expenditure,
full time equivalent staff, bed count, and total caseload by type of organization.
}

\source{
Substance Abuse and Mental Health Services Administration
}
\references{
Manderscheid, R.W. and Henderson, M.J. (2002). Mental Health, United States, 2002.

DHHS Publication No. SMA04-3938. Rockville MD USA: Substance Abuse and
Mental Health Services Administration.

Available at \url{http://mentalhealth.samhsa.gov/publications/allpubs/SMA04-3938/AppendixA.asp}
}

\seealso{
\code{\link{smho.N874}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}
}
\examples{
data(smho98)
str(smho98)
summary(smho98)
}
\keyword{datasets}
