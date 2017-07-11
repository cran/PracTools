\name{BW2stagePPSe}
\alias{BW2stagePPSe}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimated relvariance components for 2-stage sample
}
\description{
Estimate components of relvariance for a sample design where primary sampling units (PSUs) are selected with \emph{pps} and elements are selected via \emph{srs}. The input is a sample selected in this way.
}
\usage{
BW2stagePPSe(Ni, ni, X, psuID, w, m, pp)
}

\arguments{
  \item{Ni}{
 vector of number of elements in the population of each sample PSU; length is the number of PSUs in the sample.
}
  \item{ni}{
 vector of number of sample elements in each sample PSU; length is the number of PSUs in the sample. PSUs must be in the same order in \code{ni} and in \code{X}.
}
  \item{X}{
 data vector for sample elements; length is the number of elements in the sample. These must be in PSU order. PSUs must be in the same order in \code{ni} and in \code{X}.
}
  \item{psuID}{
 vector of PSU identification numbers.  This vector must be as long as \code{X}. Each element in a given PSU should have the same value in \code{psuID}.
}
  \item{w}{
 vector of full sample weights. This vector must be as long as \code{X}. Vector must be in the same order as \code{X}.
}
  \item{m}{number of sample PSUs}
  \item{pp}{
 vector of 1-draw probabilities for the PSUs. The length of this vector is the number of PSUs in the sample. Vector must be in the same order as \code{Ni} and \code{ni}.
  }
}
\details{
\code{BW2stagePPSe} computes the between and within population variance and relvariance
 components appropriate for a two-stage sample in which PSUs are selected with varying
  probabilities and with replacement. Elements within PSUs are selected by simple random sampling.
   The number of elements selected within each sample PSU can vary but must be at least two.
    The estimated components are appropriate for approximating the relvariance of the \emph{pwr}-estimator
     of a total when the same number of elements are selected within each sample PSU.
     This function can also be used if PSUs are selected by \emph{srswr} by appropriate definition of \code{pp}.
}

\value{
 List with values:
  \item{Vpsu}{estimated between PSU unit variance}
  \item{Vssu}{estimated within PSU unit variance}
  \item{B}{estimated between PSU unit relvariance}
  \item{W}{estimated within PSU unit relvariance}
  \item{k}{estimated ratio of \code{B+W} to estimated unit relvariance of the analysis variable}  
  \item{delta}{intraclass correlation estimated as \code{B/(B+W)}}
}

\references{
Cochran, W.G. (1977, pp.308-310). \emph{Sampling Techniques}. New York: John Wiley & Sons.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.4.1). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{BW2stagePPS}}, \code{\link{BW2stageSRS}}, \code{\link{BW3stagePPS}}, \code{\link{BW3stagePPSe}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
require(sampling)
require(reshape)      # has function that allows renaming variables
data(MDarea.pop)
Ni <- table(MDarea.pop$TRACT)
m <- 20
probi <- m*Ni / sum(Ni)
    # select sample of clusters
sam <- cluster(data=MDarea.pop, clustername="TRACT", size=m, method="systematic",
                pik=probi, description=TRUE)
    # extract data for the sample clusters
samclus <- getdata(MDarea.pop, sam)
samclus <- rename(samclus, c(Prob = "pi1"))


    # treat sample clusters as strata and select srswor from each
s <- strata(data = as.data.frame(samclus), stratanames = "TRACT",
            size = rep(50,m), method="srswor")
# extracts the observed data
samdat <- getdata(samclus,s)
samdat <- rename(samdat, c(Prob = "pi2"))

    # extract pop counts for PSUs in sample
pick <- names(Ni) \%in\% sort(unique(samdat$TRACT))
Ni.sam <- Ni[pick]
pp <- Ni.sam / sum(Ni)
wt <- 1/samdat$pi1/samdat$pi2

BW2stagePPSe(Ni = Ni.sam, ni = rep(50,20), X = samdat$y1,
            psuID = samdat$TRACT, w = wt,
            m = 20, pp = pp)
}
}
\keyword{methods}
\keyword{survey} 