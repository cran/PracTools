\name{deffCR}
\alias{deffCR}
\title{
Chen-Rust design effect
}
\description{
Chen-Rust design effect for an estimated mean from a stratified, clustered, two-stage samples
}
\usage{
deffCR(w, strvar=NULL, clvar=NULL, Wh=NULL, nest=FALSE, y)
}
\arguments{
  \item{w}{vector of weights for a sample; weights should be scaled for estimating population totals}
  \item{strvar}{vector of stratum identifiers; equal in length to that of \code{w}. }
  \item{clvar}{vector of cluster identifiers; equal in length to that of \code{w}. }
  \item{Wh}{vector of the proportions of elements that are in each stratum; length is number of strata.}
  \item{nest}{Are cluster IDs numbered within strata (\code{TRUE} or \code{FALSE})? If \code{TRUE}, cluster IDs can be restarted within strata, e.g., 1,2,3,1,2,3,...}
  \item{y}{vector of the sample values of an analysis variable}

}
\details{
The Chen-Rust \emph{deff} for an estimated mean accounts for stratification, clustering, and unequal weights, but does not account for the use of any auxiliary data in the estimator of a mean. The Chen-Rust \emph{deff} returned here is appropriate for stratified, two-stage sampling. Note that separate \emph{deff}'s are produced for weighting, clustering, and stratification within each stratum. These cannot be added across strata unless the stratum values of the coefficient of variation of the weights, the sample size of clusters, and the intracluster correlation of \code{y} are equal across all strata (see Chen and Rust 2017, p.117).

The weights in \code{w} should be scaled for estimating population totals. If the weights sum to the sample size or less in any stratum or overall, an error will result. This can occur if normalized weights are used.
}
\value{
  A list with components:
  \item{strata components}{Matrix with number of sample first-stage units, intracluster correlation, coefficient of variation of the weights, and \emph{deff}'s due to weighting (\code{deff.w}), clustering (\code{deff.c}), and stratification (\code{deff.s})for each stratum. When \code{strvar} or \code{clvar} are \code{NULL} appropriate subsets of these are output.}
  \item{overall deff}{Design effect for full sample accounting for weighting, clustering, and stratification}
}

\references{
Chen, S. and Rust, K. (2017). An Extension of Kish's Formula for Design Effects to Two- and Three-Stage Designs with Stratification. \emph{Journal of Survey Statistics and Methodology}, 5(2), 111-130.

Valliant, R., Dever, J., Kreuter, F. (2018, chap. 14). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{deff}}, \code{\link{deffH}}, \code{\link{deffK}}, \code{\link{deffS}}
}
\examples{
require(sampling)
require(reshape)
data(MDarea.popA)
Ni <- table(MDarea.popA$TRACT)
m <- 20
probi <- m*Ni / sum(Ni)
    # select sample of clusters
set.seed(-780087528)
sam <- sampling::cluster(data=MDarea.popA, clustername="TRACT", size=m, method="systematic",
                pik=probi, description=TRUE)
    # extract data for the sample clusters
samclus <- getdata(MDarea.popA, sam)
samclus <- rename(samclus, c("Prob" = "pi1"))
    # treat sample clusters as strata and select srswor from each
nbar <- 8
s <- sampling::strata(data = as.data.frame(samclus), stratanames = "TRACT",
            size = rep(nbar,m), method="srswor")
    # extracts the observed data
samdat <- getdata(samclus,s)
samdat <- rename(samdat, c("Prob" = "pi2"))
    # add a fake stratum ID
H <- 2
nh <- m * nbar / H
stratum <- NULL
for (h in 1:H){
    stratum <- c(stratum, rep(h,nh))
}
wt <- 1/(samdat$pi1*samdat$pi2) * runif(m*nbar)
samdat <- cbind(subset(samdat, select = -c(stratum)), stratum, wt)
deffCR(w = samdat$wt, strvar = samdat$stratum, clvar = samdat$TRACT, Wh=NULL, y=samdat$y2)
}

\keyword{methods}
\keyword{survey} 
