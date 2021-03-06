\name{BW2stageSRS}
\alias{BW2stageSRS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Relvariance components for 2-stage sample
}
\description{
Compute components of relvariance for a sample design where primary sampling units (PSUs) and elements are selected via \emph{srs}. The input is an entire sampling frame.
}
\usage{
BW2stageSRS(X, psuID)
}

\arguments{
  \item{X}{
 data vector; length is the number of elements in the population.
}
  \item{psuID}{
 vector of PSU identification numbers.  This vector must be as long as \code{X}. Each element in a given PSU should have the same value in \code{psuID}. PSUs must be in the same order as in \code{X}.
}
}
\details{
\code{BW2stageSRS} computes the between and within population relvariance components
appropriate for a two-stage sample in which PSUs are selected via \emph{srs} (either with or without replacement). Elements within PSUs are assumed to be selected by \emph{srswor}. The same number of elements is assumed to be selected within each sample PSU. The function requires that an entire frame of PSUs and elements be input.
(Use \code{\link{BW2stagePPSe}} if only a sample of PSUs and elements are available.)
}

\value{
 List with values:
  \item{B2}{between PSU unit relvariance}
  \item{W2}{within PSU unit relvariance}
  \item{unit relvar}{unit relvariance for population}
  \item{B2+W2}{\eqn{B^2 + W^2}}
  \item{k}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{delta full}{intraclass correlation estimated as \eqn{B^2/(B^2 + W^2)}}
}

\references{
Cochran, W.G. (1977, chap. 11). \emph{Sampling Techniques}. New York: John Wiley & Sons.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.2.1). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{BW2stagePPS}}, \code{\link{BW2stagePPSe}}, \code{\link{BW3stagePPS}}, \code{\link{BW3stagePPSe}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(MDarea.pop)
MDsub <- MDarea.pop[1:100000,]
    # psu's are defined by PSU variable
BW2stageSRS(abs(MDsub$Hispanic-2), psuID=MDsub$PSU)
    # psu's are defined by SSU variable
BW2stageSRS(abs(MDsub$Hispanic-2), psuID=MDsub$SSU)
}

\keyword{methods}
\keyword{survey} 