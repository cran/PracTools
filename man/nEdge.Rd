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
      alloc = NULL, Costh = NULL, nh.lowerbd = NULL
)
}
\arguments{
  \item{ci.lev}{
    desired confidence level for a 1- or 2-sided normal approximation confidence interval based on an estimated mean; must be in the interval (0,1)
}
  \item{side}{
    either \code{"two.sided"}, \code{"one.sided.upper"} or \code{"one.sided.lower"} for type of confidence interval
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
  \item{Costh}{
    vector of costs per unit in each stratum; these exclude fixed costs that do not vary with the sample size
}
  \item{nh.lowerbd}{
    Scalar value of the minimum value of the sample size allowed for any stratum  
}
}
\details{
\code{nEdge} computes the total sample size needed in either a stratified or unstratified simple random sample so that 
the coverage probability of a confidence interval is within a specified tolerance (\code{epsilon}) of a nominal confidence 
level (\code{ci.lev}). The calculation assumes that there is a single estimated mean of the variable \code{ycol} 
that is of key importance in a sample. Confidence intervals for the finite population mean are usually computed using the 
normal approximation whose accuracy depends on the underlying structure of the analytic variable and the total sample size. 

In some applications, assuring that CIs have near nominal coverage is critical. For example, for some items on business tax 
returns the US Internal Revenue Service allows sample estimates to be used but sets precision standards based on the lower 
(or upper) limit of a 1-sided CI. In some cases the computed sample size from \code{nEdge} may be too small to meet other 
goals of a sample (such as a precision requirement) and may need to be adjusted.

Using a three-term Edgeworth approximation to the distribution of the estimated overall mean in Qing & Valliant (2026), 
which extends results in Qing & Valliant (2025), \code{nEdge} computes the total sample size needed so that a CI will have coverage 
equal to the nominal value in \code{ci.lev} plus or minus 
the tolerance \code{epsilon}. The remainder term in the approximation used in \code{nEdge} is \eqn{O(n^{-3/2})}. In previous versions of 
\code{PracTools}, the remainder term in the approximation used in \code{nEdge} was \eqn{O(n^{-1})}, i.e. this version of \code{nEdge}
provides a more accurate sample size approximation. See Qing & Valliant (2026) for comparisons between the two levels of approximation.

This was also the level of approximation used in the function 
\code{nEdgeSRS} in previous versions of \code{Practools}. \code{nEdgeSRS} has been deleted from the package because \code{nEdge}
also can be used for simple random sampling. The calculation in \code{nEdge} assumes that the sampling fraction in each stratum is negligible. The total sample size 
returned by \code{nEdge} is based on the overall Edgeworth criterion; the resulting stratum sample sizes may not be large enough so 
that the normal approximation is adequate for each stratum estimator. 

When \code{dat} is a sample, the weights (\code{wts}) used in 
the estimator of the mean are assumed to be scaled for estimating population totals. They can be inverse selection probabilities, 
i.e. ones used in the \eqn{\pi}-estimator, or weights that have been adjusted to account for nonresponse or coverage errors.

The value of \code{nh.lowerbd} is usually set to be 2 or more to permit variances to be estimated.
}
\value{
  If there are no strata, the value is a list with components:
    \item{CI type}{one-sided lower, one-sided upper, or two-sided}
    \item{epsilon}{tolerance on CI coverage}
    \item{Total sample size}{numeric sample size}
    \item{allocation}{type of allocation to strata or NULL if no strata are used}
    \item{g1}{skewness of y}
    \item{g2}{excess kurtosis of y, i.e., 4th central moment divided by the square of the variance minus 3}
    \item{q1}{1st polynomial in Edgeworth expansion of the t-statistic}
    \item{q2}{2nd polynomial in Edgeworth expansion of the t-statistic}
  If the sample design has strata, the value is a list with the first four items above followed by   
    \item{Stratum values}{Data frame with columns for stratum, number of population units in each stratum, number of sample units 
            allocated to each stratum (\code{nh}), value of \code{ch} which depends on the value of \code{alloc} (see Qing & Valliant 2026), stratum skewness \code{g1h},
            stratum excess kurtosis \code{g2h}}
  If \code{alloc} = "totcost" or "totvar", the total cost of the allocation is returned and \code{Stratum values} contains a column for \code{Costh}.
}

\references{
Qing, S. and Valliant, R. (2025). Extending Cochran's Sample Size Rule to Stratified Simple Random Sampling with Applications to Audit Sampling. \emph{Journal of Official Statistics}, 41(1), 309-328.

Qing, S. and Valliant, R. (2026). Minimum Sample Size Rule for STSRS: A Three-Term Edgeworth Expansion for the Studentized Sample Mean, accepted by \emph{Journal of Survey Statistics and Methodology}. 

U.S. Internal Revenue Service (2011). 26 CFR 601.105: Examination of returns and claims for refund, credit or abatement: determination of correct tax liability. Washington DC. \url{https://www.irs.gov/pub/irs-drop/rp-11-42.pdf}
}

\author{
Richard Valliant, Siyu Qing
}
\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
}
\examples{
require(PracTools)
set.seed(1289129963)
pop <- HMT(N=10000, H=5)
    # run for full population, 2-sided CI, Neyman allocation
nEdge(ci.lev = 0.95, side = "two.sided", epsilon = 0.005, dat = pop,
      pop.sw = TRUE, wts = NULL, hcol = "strat", ycol = "y", alloc = "neyman",
      Costh = NULL)
nEdge(ci.lev = 0.95, side = "two.sided", epsilon = 0.005,dat = pop,
      pop.sw = TRUE, wts = NULL, hcol = "strat", ycol = "y", alloc = "totcost",
      Costh = c(5, 10, 15, 20, 25)
)
      
    # run for a stratified sample, 1-sided lower CI, equal allocation
require(sampling)
sam <- strata(data=pop, stratanames="strat", size=c(30, 40, 50, 60, 70), method=c("srswor"),
              description=TRUE)
samdat <- pop[sam$ID_unit,]
w = 1/sam$Prob
nEdge(ci.lev = 0.99, side = "one.sided.lower", epsilon = 0.01, dat = samdat,
      pop.sw = FALSE, wts = w, hcol = "strat", ycol = "y", alloc = "equal",
      Costh = NULL, nh.lowerbd = 2)
    
    # simple radom sample size      
nEdge(ci.lev = 0.95, side = "two.sided", epsilon = 0.005, dat = pop,
      pop.sw = TRUE, wts = NULL, hcol = NULL, ycol = "y")
}

\keyword{methods}
\keyword{survey} 
