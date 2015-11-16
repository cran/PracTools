\name{BW3stagePPS}
\alias{BW3stagePPS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relvariance components for 3-stage sample
}
\description{
Compute components of relvariance for a sample design where primary sampling units (PSUs) are selected with \emph{ppswr} and secondary sampling units (SSUs) and elements within SSUs are selected via \emph{srs}. The input is an entire sampling frame.
}
\usage{
BW3stagePPS(X, pp, psuID, ssuID)
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
  \item{ssuID}{
 vector of SSU identification numbers.  This vector must be as long as \code{X}. Each element in a given SSU should have the same value in \code{ssuID}. PSUs and SSUs must be in the same order as in \code{X}. \code{ssuID} should have the form \code{psuID}||(\code{ssuID} within PSU)}.
}
\details{
\code{BW3stagePPS} computes the between and within population relvariance components
appropriate for a three-stage sample in which PSUs are selected with varying probabilities
and with replacement. SSUs and elements within SSUs are selected by simple random sampling.
The components are appropriate for approximating the relvariance of the \emph{pwr}-estimator of a total
when the same number of SSUs are selected within each PSU, and the same number of elements are
selected within each sample SSU. The function requires that an entire sampling frame of PSUs
and elements be input. (Use \code{\link{BW2stagePPSe}} if only a sample of PSUs, SSUs, and elements is available.)
}

\value{
 List with values:
  \item{B}{between PSU unit relvariance}
  \item{W}{within PSU unit relvariance computed as if the sample were two-stage}
  \item{W2}{unit relvariance among SSU totals}
  \item{W3}{unit relvariance among elements within PSU/SSUs}
  \item{unit relvar}{unit relvariance for population}
  \item{k1}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{k2}{ratio of \eqn{W_{2}^2 + W_{3}^2} to unit relvariance}
  \item{delta1}{homogeneity measure among elements within PSUs estimated as \eqn{B^2/(B^2+W^2)}}
  \item{delta2}{homogeneity measure among elements within SSUs estimated as \eqn{W_{2}^2/(W_{2}^2 + W_{3}^2)}}
}

\references{
Hansen,M.H., Hurwitz,W.N., and Madow,W.G. (1953, chap. 9, p.211). \emph{Sample Survey Methods and Theory}, Vol.I. John Wiley & Sons.

Saerndal,C.E., Swensson, B., and Wretman, J. (1992, p.149). \emph{Model Assisted Survey Sampling}. Springer.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.2.4). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{BW2stagePPS}}, \code{\link{BW2stagePPSe}}, \code{\link{BW2stageSRS}}, \code{\link{BW3stagePPSe}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(MDarea.pop)
MDsub <- MDarea.pop[1:100000,]
M <- length(unique(MDsub$PSU))
    # srs/srs/srs design
pp.PSU <- rep(1/M,M)
BW3stagePPS(X=MDsub$y1, pp=pp.PSU, psuID=MDsub$PSU, ssuID=MDsub$SSU)
    # ppswr/srs/srs design
pp.PSU <- table(MDsub$PSU) / nrow(MDsub)
BW3stagePPS(X=MDsub$y1, pp=pp.PSU, psuID=MDsub$PSU, ssuID=MDsub$SSU)
}

\keyword{methods}
\keyword{survey} 