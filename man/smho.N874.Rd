\name{smho.N874}
\alias{smho.N874}
\docType{data}
\title{
Survey of Mental Health Organizations Data
}
\description{
Data from the 1998 Survey of Mental Health Organizations (SMHO)
}
\usage{data(smho.N874)}
\format{
  A data frame with 874 observations on the following 6 variables.
  \describe{
    \item{\code{EXPTOTAL}}{Total expenditures in 1998}
    \item{\code{BEDS}}{Total inpatient beds}
    \item{\code{SEENCNT}}{Unduplicated client/patient count seen during year}
    \item{\code{EOYCNT}}{End of year count of patients on the role}
    \item{\code{FINDIRCT}}{Hospital receives money from the state mental health agency (1=Yes; 2=No)}
    \item{\code{hosp.type}}{Hospital type
       (1 = Psychiatric;
       2 = Residential or veterans;
       3 = General;
       4 = Outpatient, partial care;
       5 = Multi-service, substance abuse)
     }
  }
}
\details{
The 1998 SMHO was conducted by the U.S. Substance Abuse and Mental Health Services Administration.
It collected data on mental health care organizations and general hospitals that provide mental health care services, with an objective to develop national and state level estimates for total expenditure,
full time equivalent staff, bed count, and total caseload by type of organization.
The population omits one extreme observation in the \code{smho98} population and
has fewer variables than \code{smho98}.  \code{smho.N874} contains observations on 874 facilities.
}

\source{
Substance Abuse and Mental Health Services Administration
}
\references{
Manderscheid, R.W. and Henderson, M.J. (2002). Mental Health, United States, 2002.
DHHS Publication No. SMA04-3938. Rockville MD USA: Substance Abuse and
Mental Health Services Administration.

}

\seealso{
\code{\link{smho98}}
}
\examples{
data(smho.N874)
str(smho.N874)
}
\keyword{datasets}
