\name{nAuditMUS}
\alias{nAuditMUS}

\title{
Sample sizes for a Monetary Unit Sample in an audit
}
\description{
Compute a sample size for an audit where the goal is to control the probability of observing only a small number of errors given an underlying error rate in the population. The sample will be selected with probabilities proportional to a measure of size (MOS). When the MOS of each record is a monetary unit, auditors refer to this as an monetary unit sampling or dollar unit sampling.
}
\usage{
nAuditMUS(MUSVar, Value.sw, CL = 0.90, Error.sw, Tol.Error, Exp.Error = 0)
}
\arguments{
  \item{MUSVar}{
    The measure of size for monetary unit sampling (MUS)
  }
  \item{Value.sw}{
    Determines whether the monetary unit sample is based on positive values, negative values, or absolute values. If \code{Value.sw} = "Positive" or "Pos", only positive values of \code{MUSVar} are used; if "Negative" or "Neg" only negative values are used; if "Absolute" or "Abs", all values of \code{MUSVar} are used after taking their absolute values.
  }
  \item{CL}{
    Probability that the sample will contain an acceptable number of errors. Auditors refer to this as 'confidence level'. The probability that the sample will contain the tolerable number of  errors or fewer is \code{1-CL}. The range of \code{CL} is 0 to 1.
  }
  \item{Error.sw}{
    Determines whether errors are based on monetary amounts or percentages, i.e., whether \code{Tol.error} is interpreted as a dollar amount (\code{Error.sw} = "Absolute" or "Amt") or as a percent (\code{Error.sw} = "Percent" or "Pct").
  }
  \item{Tol.Error}{
    The amount of error expressed as a value or a percentage that the auditor considers tolerable. If \code{Error.sw} is "Percent" or "Pct", \code{Tol.Error} is a percent between 0 and 100. If \code{Error.sw} = "Absolute" or "Amt", \code{Tol.Error} is interpreted as a dollar amount.
  }
  \item{Exp.Error}{
    The amount of error, expressed as a value or a percentage, that the auditor expects in the population. If \code{Error.sw} is "Percent" or "Pct", \code{Exp.Error} is a percent between 0 and 100. If \code{Error.sw} = "Absolute" or "Amt", \code{Exp.Error} is interpreted as a dollar amount.
  }
}
\details{
\code{nAuditMUS} computes the minimum sample size required for a given population, tolerable error rate or misstatement, and desired confidence level. If the expected error or misstatement is 0, (i.e. \code{Exp.Error} = 0), then the sample size is computed using the hypergeometric distribution where the acceptable number of deviations is 0. If the expected error is greater than 0, then sample size is computed by first calculating the maximum sample size where the number of deviations divided by the sample size is less than the expected error, then calculating the minimum sample size where the number of deviations divided by the sample size is greater than the expected error, and finally performing a straight line interpolation between these two values where the interpolated value is the specified expected error. The returned sample size calculation is the ceiling of that interpolated sample size.
}

\value{
 List object with values:
      \item{Value.Range}{Whether the MUS variable is for positive, negative, or absolute values as defined by \code{Value.sw}}
      \item{Error.Type}{Amount or Percent as defined by \code{Error.sw}}
      \item{Tol.Error.Rate}{The tolerable error expressed as a percentage of items if \code{Error.sw} = "Percent" or "Pct" or as a proportion of total monetary units otherwise}
      \item{Exp.Error.Rate}{The expected error expressed as a percentage of items if \code{Error.sw} = "Percent" or "Pct" or as a proportion of total monetary value otherwise}
      \item{Number.Records}{The population count of records in the value range based on selecting ones with positive, negative or absolute value of \code{MUSVar}}
      \item{Sample.Size}{Minimum sample size needed to meet tolerable and expected error rate requirements}
      \item{Number.HighVal}{Number of records that are high value (exceed the interval used for systematic sampling) and will be certainties in the sample}
      \item{Positive.Pop.Dollars}{The absolute value of the total dollar (or other monetary unit) amount in the population in the value range}
      \item{Conf.level}{Probability that the sample will meet MUS requirements}
      \item{Sampling.Interval}{Spacing or skip interval that would be used in a systematic probability proportional to monetary unit sampling}
}

\references{
GAO (2020). Financial Audit Manual, Volume 1, section 480.21-480.26. Washington DC; \url{https://www.gao.gov/assets/gao-18-601g.pdf}
}

\author{
George Zipf, Richard Valliant
}
\seealso{
\code{\link{nAuditAttr}}
}
\examples{
    # generate an artificial population with some negative monetary amounts
EX <- 1000
relvar <- 2
alpha <- 1/relvar
sigma <- EX * relvar
lowval <- 100                   # minimum positive X's allowed
prop.neg <- 0.05                # proportion of pop with negative values
N.neg <- floor(1000 * prop.neg)    # number of negative X's allowed
X <- rgamma(n=1000, shape=alpha, scale=sigma)
Xlow <- sort(X)[1:N.neg]
xneg <- -Xlow - lowval
xpos <- X[N.neg:length(X)]
X <- c(xneg, xpos)

nAuditMUS(X, Value.sw = "Pos", Error.sw = "Amount", Tol.Error = 180000, Exp.Error = 10000)
nAuditMUS(X, Value.sw = "Pos", Error.sw = "Pct", Tol.Error = 18, Exp.Error = 3)
nAuditMUS(X, Value.sw = "Abs", Error.sw = "Amount", Tol.Error = 180000, Exp.Error = 10000)
}

\keyword{methods}
\keyword{survey} 