\name{clusOpt2fixedPSU}
\alias{clusOpt2fixedPSU}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Optimal number of sample elements per PSU in a two-stage sample when the sample of PSUs is fixed
}
\description{
Compute the optimum number of sample elements per primary sampling unit (PSU) for a fixed set of PSUs
}
\usage{
clusOpt2fixedPSU(C1, C2, m, delta, unit.rv, k=1, CV0=NULL, tot.cost, cal.sw)
}

\arguments{
  \item{C1}{unit cost per PSU}
  \item{C2}{unit cost per element}
  \item{m}{number of sample PSU's (fixed)}
  \item{delta}{homogeneity measure}
  \item{unit.rv}{unit relvariance}
  \item{k}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{CV0}{target CV}
  \item{tot.cost}{total budget for variable costs}
  \item{cal.sw}{specify type of optimum:
              1 = find optimal \eqn{\bar{n}} for fixed total budget;
              2 = find optimal \eqn{\bar{n}} for target CV0}
}

\details{
\code{clusOpt2fixedPSU} will compute \eqn{\bar{n}_{opt}} for a two-stage sample which uses simple
random sampling at each stage or \emph{ppswr} at the first stage and \emph{srs} at the second. The PSU sample is fixed.
}
\value{
 List with values:
  \item{C1}{unit cost per PSU}
  \item{C2}{unit cost per element}
  \item{m}{number of (fixed) sample PSUs}
  \item{delta}{homogeneity measure}
  \item{unit relvar}{unit relvariance}
  \item{k}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{cost}{total budget for variable costs, \eqn{C-C_{0}} if \code{cal.sw}=1; or computed cost if \code{cal.sw}=2}
  \item{n}{optimum number of sample elements per PSU}
  \item{CV}{computed CV if \code{cal.sw}=1; or target CV if \code{cal.sw}=2}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, sect. 9.3.3).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{clusOpt2}}, \code{\link{clusOpt3}}, \code{\link{clusOpt3fixedPSU}}
}
\examples{
    # optima for a vector of budgets
clusOpt2fixedPSU(C1=500, C2=100, m=100, delta=0.05, unit.rv=2, k=1, CV0=NULL,
       tot.cost=c(100000, 500000, 10^6), cal.sw=1)
    # optima for a target CV and vector of PSU costs
clusOpt2fixedPSU(C1=c(500,1000,5000), C2=100, m=100, delta=0.05, unit.rv=2, k=1,
       CV0=0.05, tot.cost=NULL, cal.sw=2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{survey} 