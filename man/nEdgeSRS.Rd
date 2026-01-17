\name{nEdgeSRS}
\alias{nEdgeSRS}
\title{
Compute the total sample size for a simple random sample based on an Edgeworth approximation
}
\description{
Compute the total simple random sample size that is large enough to insure adequate coverage of a normal approximation confidence interval (CI) for a population mean.
}
\usage{
nEdgeSRS(ci.lev, side, epsilon = 0.005, dat, pop.sw = TRUE, wts = NULL, hcol=NULL, ycol)
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
}

\details{
\code{nEdgeSRS} computes the total sample size needed in a simple random sample so that the coverage probability of a confidence interval is within a specified tolerance (\code{epsilon}) of a nominal confidence level (\code{ci.lev}). Confidence intervals for the finite population mean are usually computed using the normal approximation whose accuracy depends on the sample size and the underlying structure of the analytic variable. In some applications, assuring that CIs have near nominal coverage is critical. For example, for some items on business tax returns the US Internal Revenue Service allows sample estimates to be used but sets precision standards based on the lower (or upper) limit of a 1-sided CI.

Using an Edgeworth approximation in Sugden, Smith, and Jones (SSJ, 2000) to the distribution of the estimated mean, \code{nEdgeSRS} computes the total sample size needed so that a CI will have coverage equal to the nominal value in \code{ci.lev} plus or minus the tolerance \code{epsilon}. Two alternatives are given: (1) a sample size from solving quadratic equation (4.4) in SSJ and (2) a modification of a rule from Cochran (1977) given in expression (4.3) of SSJ. If \code{hcol} is specified, a separate calculation is made in each stratum of the required stratum simple random sample size; thus, each stratum sample size should be adequate so that the normal approximation for each stratum estimator holds. The calculation assumes that the overall or stratum sampling fractions are negligible.

When \code{dat} is a sample, the weights (\code{wts}) used in the estimator of the mean (or total) are assumed to be scaled for estimating population totals. They can be inverse selection probabilities, i.e. ones used in the \eqn{\pi}-estimator, or weights that have been adjusted to account for nonresponse or coverage errors.

The remainder term in the approximation used in \code{nEdgeSRS} is \eqn{O(n^{-3/2})}. The function \code{nEdge} also uses a 
\eqn{O(n^{-3/2})} approximation but applies to an overall mean from a stratified simple random sample for which several different allocations 
can be specified. 
}
\value{
  List with values:
    \item{CI type}{one-sided or two-sided}
    \item{epsilon}{tolerance on CI coverage}
    \item{Total sample size}{vector of numeric sample sizes from (1) solving SSJ (2000) quadratic equation and (2) SSJ's modified Cochran rule}
    \item{g1}{overall skewness and kurtosis; returned if no strata are used}
    \item{Stratum values}{data frame with columns for stratum, number of sample units allocated to each stratum (\code{nh}) based on the SSJ quadratic rule, proportion that each quadratic-rule stratum sample is of the total sample (\code{ph}), modified Cochran sample size (\code{nh.cochran}), skewness in each stratum (\code{stratum.skewness}), and kurtosis in each stratum (\code{stratum.kurtosis}); returned if strata are used}
}

\references{
Cochran, W.G. (1977). \emph{Sampling Techniques}, 3rd edition. New York: Wiley.

Sugden, R. A., Smith, T. M. F., and Jones, R. P. (2000). Cochran's Rule for Simple Random Sampling. \emph{Journal of the Royal Statistical Society. Series B}, Vol. 62, No.4, 787-793.
doi:https://doi.org/10.1111/1467-9868.00264

U.S. Internal Revenue Service (2011). 26 CFR 601.105: Examination of returns and claims for refund, credit or abatement: determination of correct tax liability. Washington DC.
\url{https://www.irs.gov/pub/irs-drop/rp-11-42.pdf}
}

\author{
Richard Valliant
}
\seealso{
\code{\link{nCont}}, \code{\link{nEdge}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
}
\examples{
require(PracTools)
set.seed(1289129963)
#   test using HMT pop
pop <- HMT(N=10000, H=5)
    # using pop with no strata
nEdgeSRS(ci.lev=0.95, side="one.sided", dat=pop, pop.sw=TRUE, hcol=NULL, ycol="y")
    # using a sample as input
require(sampling)
sam <- strata(data=pop, stratanames="strat", size=c(30, 40, 50, 60, 70), method=c("srswor"),
              description=TRUE)
samdat <- pop[sam$ID_unit,]
w = 1/sam$Prob
nEdgeSRS(ci.lev=0.95, side="one.sided", epsilon=0.005, dat=samdat, pop.sw=FALSE, wts=w,
         hcol="strat", ycol="y")
}

\keyword{methods}
\keyword{survey} 