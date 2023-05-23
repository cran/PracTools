\name{CompMOS}
\alias{CompMOS}
\title{
   Compute a composite measure of size for domain-based two-stage sampling
}
\description{
   Compute a composite measure of size variable for domain-based sampling that accounts for desired sampling rates of domain units.
}
\usage{
   CompMOS(dsn = NULL,  psuID = NULL, n.PSU = NULL, domain = NULL, domain.req.n = NULL,
   exp.domain.rr = NULL)
}
\arguments{
   \item{dsn}{Data (sampling) frame used for Composite MOS calculations}
   \item{psuID}{PSU Cluster ID}
   \item{n.PSU}{PSU sample size}
   \item{domain}{Vector of domain variable names}
   \item{domain.req.n}{Vector of required sample size from each domain}
   \item{exp.domain.rr}{Vector of expected response rate for each domain as a percentage between 0 and 1}
}
\details{
   Two-stage samples are often selected from populations for which separate estimates are required for domains, i.e., subpopulations. Composite measures of size for selecting PSU samples with probability proportional to that size can accomplish three things:
   \enumerate{
      \item{Self-weighting samples from each of several domains}
      \item{Equal workload in each PSU, i.e., same total sample size in each PSU (across all domains)}
      \item{PSU selection probabilities that give "credit" for containing domains that are relatively rare in the population}
   }
   \code{CompMOS} computes a single composite measure of size, probability of inclusion for each PSU in the sampling frame, and within-PSU sampling rates for each domain. Additional variables regarding survey operations at the PSU domain level are also provided (see Value section below).
}
\value{
   A list with four components:
   \item{warning}{If domain sampling at the desired rate is not feasible in one of more PSUs (i.e. the domain population in the PSU is too small to meet the domain sampling requirements), a warning is included. Review \code{CompMOS.psuID} to see where the sampling is not feasible. If all PSUs pass the feasibility test, warning = "None".}
   \item{CompMOS.psuID}{A data frame containing the input psuID and domain variables, the composite measure of size, the probability of inclusion, and the PSU/domain sampling fractions,
       PSU/domain sample sizes, and a feasibility check on each PSU/domain to ensure that the PSU/domain population size is sufficient for sampling.}
   \item{CompMOS.design}{A data frame containing domain level  survey design and sample information from the input data frame and input domain requirements.}
   \item{CompMOS.Ops}{A data frame containing the number of PSUs, the sample workload, and the PSU workload.}
}
\references{
Aldworth J., Hirsch E. L., Martin P. C., Shook-Sa B. E. (2015). 2014 National Survey on Drug Use and Health sample design report. Tech. Rep. Prepared under contract no. HHSS283201300001C by RTI International, Substance Abuse and Mental Health Services Administration, \code{https://www.
samhsa.gov/data/sites/default/files/NSDUHmrbSampleDesign2014v1.pdf}

Singh, A.C. and Harter, R. (2015). Domain sample allocation within primary sampling units in designing domain-level equal probability selection methods. \emph{Survey Methodology}, 41(2), 297-314.

Valliant, R., Dever, J., Kreuter, F. (2018, sec. 10.5). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{George Zipf, Richard Valliant}
\examples{
psuID <- c(1:10)
D1 <- c(50, 50, 50, 50, 50, 70, 50, 50, 50, 50)
D2 <- c(50, 30, 90, 40, 25, 40, 80, 65, 30, 50)
dsn <- cbind.data.frame(psuID, D1, D2)
n.PSU <- 4
domain <- c("D1", "D2")
domain.req.n <- c(130, 50)
exp.domain.rr <- c(1, 1)

CompMOS(dsn = dsn, psuID  = psuID, n.PSU  = n.PSU, domain = domain,
        domain.req.n = domain.req.n, exp.domain.rr = exp.domain.rr)

  # MDarea.pop has multiple rows for each TRACT; need to summarize TRACT/Age totals
  # for input to CompMOS
data(MDarea.pop)
MDpop <- MDarea.pop[,1:8]
MDpop$AgeGrp <- cut(MDpop$Age, breaks = c(0, 12, 17, 23),
                    labels = c("Age.44.or.under", "Age.45-64", "Age.65+"))

xx <- by(MDpop$TRACT, INDICES=MDpop$AgeGrp, table)
  # All tracts do not contain every age group; merge tract/domain count tables, retaining all tracts
xx1 <- cbind(tract=rownames(xx$Age.44.or.under), as.data.frame(unname(xx$Age.44.or.under)))
colnames(xx1)[3] <- 'Age.44.or.under'
xx2 <- cbind(tract=rownames(xx$`Age.45-64`), as.data.frame(unname(xx$`Age.45-64`)))
colnames(xx2)[3] <- 'Age.45-64'
xx3 <- cbind(tract=rownames(xx$`Age.65+`), as.data.frame(unname(xx$`Age.65+`)))
colnames(xx3)[3] <- 'Age.65+'
pop <- merge(xx1,xx2,by='tract', all=TRUE)
pop <- merge(pop,xx3,by='tract', all=TRUE)
pop <- pop[, -c(2,4,6)]
  # recode counts for missing tract/age-groups to 0
pop[is.na(pop)] <- 0

   # Note that one tract cannot be sampled at the desired rate for the 'Age.65+' domain
MDmos <- CompMOS(dsn  = pop, psuID  = pop$tract, n.PSU  = 15,
              domain = c("Age.44.or.under", "Age.45-64", "Age.65+"),
              exp.domain.rr = c(0.60, 0.70, 0.85),
              domain.req.n  = c(100, 100, 100))
}
\keyword{methods}
\keyword{survey}
