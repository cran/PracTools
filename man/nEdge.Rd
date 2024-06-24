\name{nEdge}
\alias{nEdge}
\title{
Compute the total sample size for a stratified, simple random sample based on an Edgeworth approximation
}
\description{
Compute the total stratified, simple random sample size for various allocations that is large enough to insure adequate coverage of a normal approximation confidence interval (CI) for a population mean.
}
\usage{
nEdge(ci.lev, side, epsilon = 0.005, dat, pop.sw = TRUE, wts = NULL, hcol=NULL, ycol,
      alloc = NULL, Ch = NULL)
}
\arguments{
  \item{ci.lev}{
    desired confidence level for a 1- or 2-sided normal approximation confidence interval based on an estimated mean; must be in the interval (0,1)
}
  \item{side}{
    either \code{"two.sided"} or \code{"one.sided"} for type of confidence interval
}
  \item{epsilon}{
    tolerance on coverage probability; the sample should be large enough that CI coverage is within \eqn{\pm} \code{epsilon} of \code{ci.lev}; must be in the interval (0,1)
}
  \item{dat}{
    either a population or sample data frame
}
  \item{pop.sw}{
    TRUE if \code{dat} is for a full population; FALSE if \code{dat} is for a sample
}
  \item{wts}{
    vector of weights if \code{dat} is a sample; if \code{dat} is for a population, \code{wts = NULL}
}
  \item{hcol}{
    column of \code{dat} that contains the stratum ID; strata can be character or numeric
}
  \item{ycol}{
    column of \code{dat} that contains the analysis variable; must be numeric
}
  \item{alloc}{
    allocation to the strata; must be one of \code{prop}, \code{equal}, \code{neyman}, \code{totcost}, \code{totvar}, or \code{NULL}
}
  \item{Ch}{
    vector of costs per unit in each stratum; these exclude fixed costs that do not vary with the sample size
}
}
\details{
\code{nEdge} computes the total sample size needed in either a stratified or unstratified simple random sample so that the coverage probability of a confidence interval is within a specified tolerance (\code{epsilon}) of a nominal confidence level (\code{ci.lev}). The calculation assumes that there is a single estimated mean or total of the variable \code{ycol} that is of key importance in a sample. Confidence intervals for the finite population mean are usually computed using the normal approximation whose accuracy depends on the underlying structure of the analytic variable and the total sample size. In some applications, assuring that CIs have near nominal coverage is critical. For example, for some items on business tax returns the US Internal Revenue Service allows sample estimates to be used but sets precision standards based on the lower (or upper) limit of a 1-sided CI.

Using an Edgeworth approximation to the distribution of the estimated overall mean in Qing & Valliant (2024), \code{nEdge} computes the total sample size needed so that a CI will have coverage equal to the nominal value in \code{ci.lev} plus or minus the tolerance \code{epsilon}. The calculation assumes that the sampling fraction in each stratum is negligible. The total sample size returned by \code{nEdge} is based on the overall Edgeworth criterion; the resulting stratum sample sizes may not be large enough so that the normal approximation is adequate for each stratum estimator. When \code{dat} is a sample, the weights (\code{wts}) used in the estimator of the mean (or total) are assumed to be scaled for estimating population totals. They can be inverse selection probabilities, i.e. ones used in the \eqn{\pi}-estimator, or weights that have been adjusted to account for nonresponse or coverage errors.

The remainder term in the approximation used in \code{nEdge} is \eqn{O(n^{-1/2})}. In contrast, the function \code{nEdgeSRS} uses a \eqn{O(n^{-1})} approximation but applies only to simple random sampling.
}
\value{
  List with values:
    \item{CI type}{one-sided or two-sided}
    \item{epsilon}{tolerance on CI coverage}
    \item{Total sample size}{numeric sample size}
    \item{allocation}{type of allocation to strata or NULL if no strata are used}
    \item{Stratum values}{Data frame with columns for stratum, number of sample units allocated to each stratum (\code{nh}), proportion of sample allocated to each stratum (\code{ph}), and skewness in each stratum (\code{g1h}); if no strata are used, only \code{g1}, the overall skewness is returned}
}

\references{
Qing, S. and Valliant, R. (2024). Extending Cochran's Sample Size Rule to Stratified Simple Random Sampling with Applications to Audit Sampling. \emph{Journal of Official Statistics}, accepted.

U.S. Internal Revenue Service (2011). 26 CFR 601.105: Examination of returns and claims for refund, credit or abatement: determination of correct tax liability. Washington DC. \url{https://www.irs.gov/pub/irs-drop/rp-11-42.pdf}
}

\author{
Richard Valliant, Siyu Qing
}
\seealso{
\code{\link{nCont}}, \code{\link{nEdgeSRS}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
}
\examples{
require(PracTools)
set.seed(1289129963)
pop <- HMT(N=10000, H=5)
    # run for full population
nEdge(ci.lev=0.95, side="one.sided", dat=pop, pop.sw=TRUE, wts=NULL, hcol="strat", ycol="y",
       alloc="neyman")
    # run for a stratified sample
require(sampling)
sam <- strata(data=pop, stratanames="strat", size=c(30, 40, 50, 60, 70), method=c("srswor"),
              description=TRUE)
samdat <- pop[sam$ID_unit,]
w = 1/sam$Prob
nEdge(ci.lev=0.95, side="two.sided", epsilon=0.02, dat=samdat, pop.sw=FALSE, wts=w,
       hcol="strat", ycol="y", alloc="equal")
}

\keyword{methods}
\keyword{survey} 