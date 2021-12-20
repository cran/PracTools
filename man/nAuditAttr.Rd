\name{nAuditAttr}
\alias{nAuditAttr}

\title{
Sample sizes for an attribute sample in an audit
}
\description{
Compute a sample size for an audit where the goal is to control the probability of observing only a small number of errors given an underlying error rate in the population. Auditors refer to this as an attribute sample.
}
\usage{
nAuditAttr(TolRate=0.05, AccDev, CL, N=5000)
}
\arguments{
  \item{TolRate}{
    Proportion of units in the population with an attribute, e.g., errors in an audit. Auditors term this the 'tolerable rate of deviation' in the population to be tested.
  }
  \item{AccDev}{
    Acceptable deviation, which is the number of units with the attribute (i.e., the number of errors) that would be acceptable in the sample. The largest proportion of errors that would be deemed to be acceptable in an audit would be \code{AccDev/N}.
  }
  \item{CL}{
    Probability that the sample will contain an acceptable number of errors. Auditors refer to this as 'confidence level'. The probability that the sample will contain \code{AccDev} errors or fewer is \code{1-CL}.
  }
  \item{N}{
    Size of the population of records to be audited.
  }
}
\details{
\code{nAuditAttr} computes the minimum sample size required so that the probability, \code{1-CL} of detecting less than or equal to a specified number of errors in the sample, is controlled. The sample is assumed to be selected with equal probabilities. \code{AccDev} is the largest number of errors in the sample that will be considered as meeting the audit standards. \code{TolRate} is the underlying population error rate, which is typically set to be larger than \code{AccDev/N}. The sample size is computed in two ways: (1) using the hypergeometric distribution, which accounts for the size of the population and (2) with the binomial distribution, which will be appropriate if the population is very large. When \code{N} is large and the sampling fraction is small, both sample sizes will be approximately the same.
}

\value{
 List object with values:
      \item{"Pop.Size"}{population size}
      \item{"Tol.Dev.Rate"}{proportion of records with errors in population}
      \item{"Acceptable.Errors"}{largest number of errors, found in the sample, that will meet audit standards}
      \item{"Sample.Size.Hypergeometric"}{minimum sample size needed to detect \code{AccDev} errors in the sample computed via the hypergeometric distribution}
      \item{"Sample.Size.Binomial"}{minimum sample size needed to detect \code{AccDev} errors in the sample computed via the binomial distribution}
}

\references{
GAO (2020). Financial Audit Manual, Volume 1, section 450.08. Washington DC; \url{https://www.gao.gov/assets/gao-18-601g.pdf}

Stewart, Trevor R. (2012). \emph{Technical Notes on the AICPA Audit Guide: Audit Sampling}. American Institute of Certified Public Accountants, Inc. New York, NY 10036-8775; \url{https://us.aicpa.org/content/dam/aicpa/publications/accountingauditing/keytopics/downloadabledocuments/sampling_guide_technical_notes.pdf}
}

\author{
George Zipf, Richard Valliant
}
\examples{
# Examples from the US GAO Financial Audit Manual (2020), Figure 450.1, Table 1
nAuditAttr(AccDev = 0, CL = .90)
nAuditAttr(AccDev = 1, CL = .90)
nAuditAttr(AccDev = 2, CL = .90)
nAuditAttr(AccDev = 3, CL = .90)
nAuditAttr(AccDev = 4, CL = .90)
}

\keyword{methods}
\keyword{survey} 