\name{nDomain}
\alias{nDomain}

\title{
Compute a simple random sample size for an estimated mean or total for a domain
}
\description{
Compute a simple random sample size using either a target coefficient of variation, \eqn{CV_{0}(d)}, or target variance, \eqn{V_{0}(d)}, for an estimated mean or total for a domain.
}
\usage{
nDomain(CV0d=NULL, V0d=NULL, S2d=NULL, ybarUd=NULL, N=Inf, CVpopd=NULL, Pd, est.type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CV0d}{
 target value of coefficient of variation of estimated domain mean or total
    }
  \item{V0d}{
target value of variance of estimated domain mean or total
    }
  \item{S2d}{
 unit (population) variance for domain units
    }
  \item{ybarUd}{
 population mean of target variable for domain units
    }
  \item{N}{
number of units in full finite population (not just the domain population)
    }
  \item{CVpopd}{
  unit (population) coefficient of variation for domain units
    }
  \item{Pd}{
  proportion of units in the population that are in the domain
    }
  \item{est.type}{
  type of estimate; allowable values are \code{"mean"} or  \code{"total"}
    }
}
\details{
If \code{CV0d} is the desired target, then the unit CV, \code{CVpopd},
or the domain population mean and variance, \code{ybarUd} and \code{S2d}, must also be provided.
If \code{V0d} is the constrained value, then \code{ybarUd} must be also be included in the function call. \code{CV0d}
will then be computed as \code{sqrt(V0d)/ybarUd}.
}
\value{
  numeric sample size
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, sec. 3.5.2). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nDomain(CV0d=0.05, N=Inf, CVpopd=1, Pd=0.5, est.type="total")
nDomain(CV0d=0.05, N=Inf, CVpopd=1, Pd=0.5, est.type="mean")
nDomain(V0d=50, ybarUd=50, S2d=100, N=Inf, Pd=0.5, est.type="total")
nDomain(CV0d=0.05, ybarUd=50, S2d=100, N=Inf, Pd=0.5, est.type="total")
nDomain(CV0d=0.05, ybarUd=50, S2d=100, N=Inf, Pd=0.5, est.type="mean")
}

\keyword{methods}
\keyword{survey}
