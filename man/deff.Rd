\name{deff}
\alias{deff}
\title{
Design effects of various types
}
\description{
Compute the Kish, Henry, Spencer, or Chen-Rust design effects.
}
\usage{
deff(w, x=NULL, y=NULL, p=NULL, strvar=NULL, clvar=NULL, Wh=NULL, nest=FALSE, type)
}
\arguments{
  \item{w}{vector of weights for a sample}
  \item{x}{matrix of covariates used to construct a GREG estimator of the total of \eqn{y}.  This matrix does not include the intercept. Used only for Henry \emph{deff}.}
  \item{y}{vector of the sample values of an analysis variable}
  \item{p}{vector of 1-draw selection probabilities, i.e., the probability that each unit would be selected in a sample of size 1. Used only for Spencer \emph{deff}.}
  \item{strvar}{vector of stratum identifiers; equal in length to that of \code{w}. Used only for Chen-Rust \emph{deff}. }
  \item{clvar}{vector of cluster identifiers; equal in length to that of \code{w}. Used only for Chen-Rust \emph{deff}. }
  \item{Wh}{vector of the proportions of elements that are in each stratum; length is number of strata. Used only for Chen-Rust \emph{deff}.}
  \item{nest}{Are cluster IDs numbered within strata (\code{TRUE} or \code{FALSE})? If \code{TRUE}, cluster IDs can be restarted within strata, e.g., 1,2,3,1,2,3,...}
  \item{type}{type of allocation; must be one of \code{"kish"}, \code{"henry"}, \code{"spencer"}, \code{"cr"}}
}
\details{
\code{deff} calls one of \code{deffK}, \code{deffH}, \code{deffS}, or \code{deffCR} depending on the value of the \code{type} parameter. The Kish design effect is the ratio of the variance of an estimated mean in stratified simple random sampling without replacement (\emph{stsrswor}) to the variance of the estimated mean in \emph{srswor}, assuming that all stratum unit variances are equal. In that case, proportional allocation with equal weighting is optimal. deffK equals 1 + relvar(w) where relvar is relvariance of the vector of survey weights. This measure is not appropriate in samples where unequal weighting is more efficient than equal weighting.

The Henry design effect is the ratio of the variance of the general regression (GREG) estimator of a total of \eqn{y} to the variance of the estimated total in \emph{srswr}. Calculations for the Henry \emph{deff} are done as if the sample is selected in a single-stage and with replacement.  Varying selection probabilities can be used. The model for the GREG is assumed to be \eqn{y = \alpha + \beta x + \epsilon}, i.e., the model has an intercept.

The Spencer design effect is the ratio of the variance of the \emph{pwr}-estimator of the total of \emph{y}, assuming that a single-stage sample is selected with replacement, to the variance of the total estimated in \emph{srswr}. Varying selection probabilities can be used.

The Chen-Rust \emph{deff} accounts for stratification, clustering, and unequal weights, but does not account for the use of any auxiliary data in the estimator of a mean. The Chen-Rust \emph{deff} returned here is appropriate for stratified, two-stage sampling.}
\value{
  Numeric design effect for types \code{kish}, \code{henry}, \code{spencer}. For type \code{cr} a list with components:
  \item{strata components}{Matrix with \emph{deff}'s due to weighting, clustering, and stratification for each stratum}
  \item{overall deff}{Design effect for full sample accounting for weighting, clustering, and stratification}
}

\references{
Chen, S. and Rust, K. (2017). An Extension of Kish's Formula for Design Effects to Two- and Three-Stage Designs with Stratification. \emph{Journal of Survey Statistics and Methodology}, 5(2), 111-130.

Henry, K.A., and Valliant, R. (2015). A Design Effect Measure for Calibration Weighting in Single-stage Samples. \emph{Survey Methodology}, 41, 315-331.

Kish, L. (1965). \emph{Survey Sampling}. New York: John Wiley & Sons.

Kish, L. (1992). Weighting for unequal Pi. \emph{Journal of Official Statistics}, 8, 183-200.

Park, I., and Lee, H. (2004). Design Effects for the Weighted Mean and Total Estimators under Complex Survey Sampling. \emph{Survey Methodology}, 30, 183-193.

Spencer, B. D. (2000). An Approximate Design Effect for Unequal Weighting When Measurements May Correlate With Selection Probabilities. Survey Methodology, 26, 137-138.

Valliant, R., Dever, J., Kreuter, F. (2018, chap. 14). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{deffK}}, \code{\link{deffH}}, \code{\link{deffS}}, \code{\link{deffCR}}
}
\examples{
require(reshape)      # has function that allows renaming variables
require(sampling)

set.seed(-500398777)
    # generate population using HMT function
pop.dat <- as.data.frame(HMT())
mos <- pop.dat$x
pop.dat$prbs.1d <- mos / sum(mos)
    # select pps sample
n <- 80
pk <- pop.dat$prbs.1d
sam <- UPrandomsystematic(pk)
sam <- sam==1

sam.dat <- pop.dat[sam, ]
dsgn.wts <- 1/pk[sam]
deff(w=dsgn.wts, type="kish")
deff(w=dsgn.wts, y=sam.dat$y, p=sam.dat$prbs.1d, type="spencer")
deff(w=dsgn.wts, x=sam.dat$x, y=sam.dat$y, type="henry")


data(MDarea.popA)
Ni <- table(MDarea.popA$TRACT)
m <- 10
probi <- m*Ni / sum(Ni)
    # select sample of clusters
set.seed(-780087528)
sam <- cluster(data=MDarea.popA, clustername="TRACT", size=m, method="systematic",
                pik=probi, description=TRUE)
    # extract data for the sample clusters
samclus <- getdata(MDarea.popA, sam)
samclus <- rename(samclus, c("Prob" = "pi1"))
    # treat sample clusters as strata and select srswor from each
nbar <- 4
s <- strata(data = as.data.frame(samclus), stratanames = "TRACT",
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
samdat <- cbind(subset(samdat, select = -c(Stratum)), stratum, wt)
deff(w = samdat$wt, y=samdat$y2, strvar = samdat$stratum, clvar = samdat$TRACT, Wh=NULL, type="cr")
}

\keyword{methods}
\keyword{survey} 