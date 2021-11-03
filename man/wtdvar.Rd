\name{wtdvar}
\alias{wtdvar}
\title{
Compute weighted variance
}
\description{
Compute an estimate of a population unit variance from a complex sample with survey weights.
}
\usage{
wtdvar(x,w,na.rm=TRUE)
}
\arguments{
  \item{x}{
 data vector
}
  \item{w}{
 vector of survey weights; must be same length as \code{x}
}
  \item{na.rm}{
 remove missing values (TRUE or FALSE)
}
}
\details{
 \code{wtdvar} is also used by \code{\link{BW3stagePPSe}} in estimating relvariance components.
}

\value{
  numeric estimate of population unit variance
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 9). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. Springer.
}

\author{
    Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\examples{
x <- c(1:3)
wts <- c(4, 6, 8)
wtdvar(x=x, w=wts)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{survey}
