\name{wtd.moments}
\alias{wtd.moments}
\title{
Compute moments of a variable from either a population or sample
}
\description{
Compute the 2nd, 3rd, 4th moments, skewness, and kurtosis of a variable from either population or sample input
}
\usage{
wtd.moments(y, w=NULL, pop.sw=TRUE)
}
\arguments{
  \item{y}{variable to be analyzed}
  \item{w}{vector of weights if the input is a sample}
  \item{pop.sw}{is the input for a population (\code{pop.sw=TRUE}) or for a sample (\code{pop.sw=FALSE})?}
}

\details{
The \eqn{r^{th}} population moment is defined as \eqn{m_r = (1/N) \sum_{k \in U} (y_k - \bar{y}_U)^r} where \emph{U} is the set of population units, 
\emph{N} is the population size, and \eqn{\bar{y}_U} is the population mean. 
When the input is for the whole population, \code{wtd.moments} evaluates this directly for \eqn{r=2, 3, 4}. 
When the input is for a sample, the \eqn{r^{th}} moment is estimated as \eqn{\hat{m}_r = (K/\hat{N}) \sum_{k \in s} ( w_k (y_k - \hat{\bar{y}}_U)^r ), r=2, 3, 4} 
where \eqn{s} is the set of sample units, \eqn{w_k} is the weight for sample unit \eqn{k}, \eqn{\hat{N} = \sum_s w_k}, and 
\eqn{\hat{\bar{y}}_U = \sum_{k \in s} w_k y_k / \hat{N}}. 
When \eqn{r=2}, \eqn{K=n/(n-1)} so that the estimator equals the unbiased variance estimator if the sample is a simple random sample; 
if \eqn{r=3,4}, then \eqn{K=1}. The function also computes or estimates the population skewness, defined as \eqn{m_3/m_2^{3/2}} 
and the population kurtosis, \eqn{m_4/m_2^2}.

The weights should be scaled for estimating population totals. The sample can be obtained from any complex design. The estimators produced by \code{wtd.moments} are consistent for general sample designs but not unbiased (except for \eqn{r=2} and SRS).
}
\value{
 Vector with values:
  \item{m2}{2nd moment}
  \item{m3}{3rd moment}
  \item{m4}{4th moment}
  \item{skewness}{skewness}
  \item{kurtosis}{kurtosis}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, sect. 3.4).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{wtdvar}}
}
\examples{
require(PracTools)
wtd.moments(y = hospital$y, w = NULL)
require(sampling)
sam <- strata(data = labor, stratanames = "h", size = c(30, 20, 10), method = c("srswor"),
              description=TRUE)
samdat <- labor[sam$ID_unit,]
wtd.moments(y = samdat$WklyWage, w = 1/sam$Prob, pop.sw=FALSE)
}

\keyword{methods}
\keyword{survey} 
