\name{BW2stagePPS}
\alias{BW2stagePPS}

\title{
Relvariance components for 2-stage sample
}
\description{
Compute components of relvariance for a sample design where primary sampling units (PSUs) are selected with probability proportional to size (\emph{pps}) and elements are selected via simple random sampling (\emph{srs}). The input is an entire sampling frame.
}
\usage{
BW2stagePPS(X, pp, psuID)
}

\arguments{
  \item{X}{
 data vector; length is the number of elements in the population.
}
  \item{pp}{
 vector of one-draw probabilities for the PSUs; length is number of PSUs in population.
}
  \item{psuID}{
 vector of PSU identification numbers.  This vector must be as long as \code{X}. Each element in a given PSU should have the same value in \code{psuID}. PSUs must be in the same order as in \code{X}.
}
}
\details{
\code{BW2stagePPS} computes the between and within population relvariance components
 appropriate for a two-stage sample in which PSUs are selected with varying probabilities
  and with replacement. Elements within PSUs are selected by simple random sampling.
   The components are appropriate for approximating the relvariance of the probability-with-replacement (\emph{pwr})-estimator of a total when the same number of elements are selected within each sample PSU.
   The function requires that an entire frame of PSUs and elements be input.
   (Use \code{\link{BW2stagePPSe}} if only a sample of PSUs and elements are available.)
}

\value{
 List object with values:
  \item{B2}{between PSU unit relvariance}
  \item{W2}{within PSU unit relvariance}
  \item{unit relvar}{unit relvariance for population}
  \item{B2+W2}{sum of between and within relvariance estimates}
  \item{k}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{delta}{measure of homogeneity with PSUs estimated as \eqn{B^2/(B^2 + W^2)}}
}

\references{
Cochran, W.G. (1977, pp.308-310). \emph{Sampling Techniques}. New York: John Wiley & Sons.

S{\"a}aerndal,C.E., Swensson, B., and Wretman, J. (1992). \emph{Model Assisted Survey Sampling}. New York: Springer.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.3.2). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.

}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{BW2stagePPSe}}, \code{\link{BW2stageSRS}}, \code{\link{BW3stagePPS}}, \code{\link{BW3stagePPSe}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(MDarea.pop)
MDsub <- MDarea.pop[1:100000,]
    # Use PSU and SSU variables to define psu's
pp.PSU <- table(MDsub$PSU) / nrow(MDsub)
pp.SSU <- table(MDsub$SSU) / nrow(MDsub)
    # components with psu's defined by the PSU variable
BW2stagePPS(MDsub$y1, pp=pp.PSU, psuID=MDsub$PSU)
    # components with psu's defined by the SSU variable
BW2stagePPS(MDsub$y1, pp=pp.SSU, psuID=MDsub$SSU)

    # Use census tracts and block groups to define psu's
trtBG <- 10*MDsub$TRACT + MDsub$BLKGROUP
pp.trt <- table(MDsub$TRACT) / nrow(MDsub)
pp.BG <- table(trtBG) / nrow(MDsub)
    # components with psu's defined by tracts
BW2stagePPS(MDsub$ins.cov, pp=pp.trt, psuID=MDsub$TRACT)
    # components with psu's defined by block groups
BW2stagePPS(MDsub$ins.cov, pp=pp.BG, psuID=trtBG)
}

\keyword{methods}
\keyword{survey} 