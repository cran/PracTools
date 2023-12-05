\name{BW3stagePPSe}
\alias{BW3stagePPSe}

\title{
Estimated relvariance components for 3-stage sample
}
\description{
Estimate components of relvariance for a sample design where primary sampling units (PSUs) are selected with probability proportional to size with replacement (\emph{ppswr}) and secondary sampling units (SSUs) and elements within SSUs are selected via simple random sampling (\emph{srs}). The input is a sample selected in this way.
}
\usage{
BW3stagePPSe(dat, v, Ni, Qi, Qij, m, lonely.SSU = "mean", lonely.TSU = "mean")
}

\arguments{
  \item{dat}{
 data frame for sample elements with PSU and SSU identifiers, weights, and analysis variable(s).  The data frame should be sorted in hierarchical order: by PSU and SSU within PSU.
    Required names for columns:
        \code{psuID} = PSU identifier;
        \code{ssuID} = SSU identifier. These must be unique, i.e., numbering should not restart within each PSU. Setting \code{ssuID} = \code{psuID}||(\code{ssuID} within PSU) is a method of doing this.
        \code{w1i}  = vector of weights for PSUs;
        \code{w2ij} = vector of weights for SSUs (PSU weight*SSU weight within PSU);
        \code{w} = full sample weight
}
  \item{v}{
 Name or number of column in data frame \code{dat} with variable to be analyzed.
}
  \item{Ni}{
 \code{m}-vector of number of SSUs in the population in the sample PSUs; \code{m} is number of sample PSUs.
}
  \item{Qi}{
 \code{m}-vector of number of elements in the population in the sample PSUs
}
  \item{Qij}{
 vector of numbers of elements in the population in the sample SSUs
}
  \item{m}{number of sample PSUs}
  \item{lonely.SSU}{
 indicator for how singleton SSUs should be handled when computing the within PSU unit relvariance. Allowable values are \code{"mean"} and \code{"zero"}.
}
 \item{lonely.TSU}{
 indicator for how singleton third-stage units (TSUs) should be handled when computing the within SSU unit relvariance. Allowable values are \code{"mean"} and \code{"zero"}.
}
}
\details{
\code{BW3stagePPSe} computes the between and within population relvariance components appropriate
for a three-stage sample in which PSUs are selected with varying probabilities and with replacement.
SSUs and elements within SSUs are selected by simple random sampling.
The estimated components are appropriate for approximating the relvariance of the
\emph{pwr}-estimator of a total when the same number of SSUs are selected within each PSU,
and the same number of elements are selected within each sample SSU.

If a PSU contains multiple SSUs, some of which have missing data, or contains only one SSU, a value is imputed. If \code{lonely.SSU = "mean"}, the mean of the non-missing PSU contributions is imputed. If \code{lonely.SSU = "zero"}, a 0 is imputed. The former would be appropriate if a PSU contains multiple SSUs but one or more of them has missing data in which case R will normally calculate an NA. The latter would be appropriate if the PSU contains only one SSU which would be selected with certainty in any sample. \code{lonely.TSU} has a similar purpose for third-stage units.
}

\value{
 List with values:
  \item{Vpsu}{estimated between PSU unit variance}
  \item{Vssu}{estimated second-stage unit variance among SSU totals}
  \item{Vtsu}{estimated third-stage unit variance}
  \item{B}{estimated between PSU unit relvariance}
  \item{W}{estimated within PSU unit relvariance computed as if the sample were two-stage}
  \item{k1}{estimated ratio of \code{B+W} to estimated unit relvariance of the analysis variable}
  \item{W2}{estimated unit relvariance among SSU totals}
  \item{W3}{estimated third-stage unit relvariance among elements within PSU/SSUs}
  \item{k2}{estimated ratio of \code{W2+W3} to estimated unit relvariance of the analysis variable}
  \item{delta1}{homogeneity measure among elements within PSUs estimated as \eqn{B^2/(B^2+W^2)}}
  \item{delta2}{homogeneity measure among elements within SSUs estimated as \eqn{W_{2}^2/(W_{2}^2 + W_{3}^2)}}
}

\references{
Hansen, M.H., Hurwitz, W.N., and Madow, W.G. (1953, chap. 9, sect. 10). \emph{Sample Survey Methods and Theory}, Vol.II. New York: John Wiley & Sons.

Valliant, R., Dever, J., Kreuter, F. (2018, sect. 9.4.2). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{BW2stagePPS}}, \code{\link{BW2stagePPSe}}, \code{\link{BW2stageSRS}}, \code{\link{BW3stagePPS}}
}
\examples{
\donttest{
    # select 3-stage sample from Maryland population
set.seed(-780087528)
data(MDarea.popA)
MDpop <- MDarea.popA
require(sampling)
require(plyr)      # has function that allows renaming variables
    # make counts of SSUs and elements per PSU
xx <- do.call("rbind",list(by(1:nrow(MDpop),MDpop$SSU,head,1)))
pop.tmp <- MDpop[xx,]
Ni <- table(pop.tmp$PSU)
Qi <- table(MDarea.popA$PSU)
Qij <- table(MDpop$SSU)
m <- 30         # no. of PSUs to select
probi <- m*Qi / sum(Qi)
    # select sample of clusters
sam <- cluster(data=MDpop, clustername="PSU", size=m, method="systematic",
               pik=probi, description=TRUE)
    # extract data for the sample clusters
samclus <- getdata(MDpop, sam)
samclus <- rename(samclus, c("Prob" = "p1i"))
samclus <- samclus[order(samclus$PSU),]
    # treat sample clusters as strata and select srswor of block groups from each
    # identify psu IDs for 1st instance of each ssuID
xx <- do.call("rbind",list(by(1:nrow(samclus),samclus$SSU,head,1)))
SSUs <- cbind(PSU=samclus$PSU[xx], SSU=samclus$SSU[xx])
    # select 2 SSUs per tract
n <- 2
s <- strata(data = as.data.frame(SSUs), stratanames = "PSU",
            size = rep(n,m), method="srswor")
s <- rename(s, c("Prob" = "p2i"))
    # extract the SSU data
    # s contains selection probs of SSUs, need to get those onto data file
SSUsam <- SSUs[s$ID_unit, ]
SSUsam <- cbind(SSUsam, s[, 2:3])
    # identify rows in PSU sample that correspond to sample SSUs
tmp <- samclus$SSU \%in\% SSUsam$SSU
SSUdat <- samclus[tmp,]
SSUdat <- merge(SSUdat, SSUsam[, c("p2i","SSU")], by="SSU")
    # select srswor from each sample SSU
n.SSU <- m*n
s <- strata(data = as.data.frame(SSUdat), stratanames = "SSU",
            size = rep(50,n.SSU), method="srswor")
s <- rename(s, c("Prob" = "p3i"))
samclus <- getdata(SSUdat, s)
del <- (1:ncol(samclus))[dimnames(samclus)[[2]] \%in\% c("ID_unit","Stratum")]
samclus <- samclus[, -del]
    # extract pop counts for PSUs in sample
pick <- names(Qi) \%in\% sort(unique(samclus$PSU))
Qi.sam <- Qi[pick]
    # extract pop counts of SSUs for PSUs in sample
pick <- names(Ni) \%in\% sort(unique(samclus$PSU))
Ni.sam <- Ni[pick]
    # extract pop counts for SSUs in sample
pick <- names(Qij) \%in\% sort(unique(samclus$SSU))
Qij.sam <- Qij[pick]
    # compute full sample weight and wts for PSUs and SSUs
wt <- 1 / samclus$p1i / samclus$p2i / samclus$p3i
w1i <- 1 / samclus$p1i
w2ij <- 1 / samclus$p1i / samclus$p2i
samdat <- data.frame(psuID = samclus$PSU, ssuID = samclus$SSU,
                     w1i = w1i, w2ij = w2ij, w = wt,
                     samclus[, c("y1","y2","y3","ins.cov", "hosp.stay")])
BW3stagePPSe(dat=samdat, v="y1", Ni=Ni.sam, Qi=Qi.sam, Qij=Qij.sam, m,
    lonely.SSU = "mean", lonely.TSU = "mean")
}
}

\keyword{methods}
\keyword{survey} 