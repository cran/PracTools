\name{clusOpt2}
\alias{clusOpt2}
\title{
Compute optimal sample sizes for a two-stage sample
}
\description{
Compute the sample sizes that minimize the variance of the \emph{pwr}-estimator of a total in a two-stage sample.
}
\usage{
clusOpt2(C1, C2, delta, unit.rv, k=1, CV0=NULL, tot.cost=NULL, cal.sw)
}
\arguments{
  \item{C1}{unit cost per primary sampling unit (PSU)}
  \item{C2}{unit cost per element}
  \item{delta}{homogeneity measure \eqn{\delta}}
  \item{unit.rv}{unit relvariance}
  \item{k}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{CV0}{target CV}
  \item{tot.cost}{total budget for variable costs}
  \item{cal.sw}{specify type of optimum:
              1 = find optimal \code{m.opt} for fixed total budget;
              2 = find optimal \code{m.opt} for target \code{CV0}}
}

\details{
\code{clusOpt2} will compute \eqn{m_{opt}} and  \eqn{\bar{n}_{opt}} for a two-stage sample
 which uses simple random sampling at each stage or \emph{ppswr} at the first stage and \emph{srs} at the second.
}
\value{
 List with values:
  \item{C1}{unit cost per PSU}
  \item{C2}{unit cost per element}
  \item{delta}{homogeneity measure}
  \item{unit relvar}{unit relvariance}
  \item{k}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{cost}{total budget for variable costs, \eqn{C-C_{0}} if \code{cal.sw}=1; or computed cost if \code{cal.sw}=2}
  \item{m.opt}{optimum number of sample PSUs}
  \item{n.opt}{optimum number of sample elements per PSU}
  \item{CV}{computed CV if \code{cal.sw}=1; or target CV if \code{cal.sw}=2}
}

\references{
Hansen,M.H., Hurwitz,W.N., and Madow,W.G. (1953, chap. 6, sect. 16).  \emph{Sample Survey Methods and Theory}, Vol.I. John Wiley & Sons.

Valliant, R., Dever, J., Kreuter, F. (2018, sect. 9.3.1).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{clusOpt2fixedPSU}}, \code{\link{clusOpt3}}, \code{\link{clusOpt3fixedPSU}}

%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
    # optimum for a fixed total budget
clusOpt2(C1=750, C2=100, delta=0.05, unit.rv=1, k=1, tot.cost=100000, cal.sw=1)
clusOpt2(C1=750, C2=100, delta=seq(0.05,0.25,0.05), unit.rv=1, k=1, tot.cost=100000, cal.sw=1)
    # optimum for a target CV
clusOpt2(C1=750, C2=100, delta=0.01, unit.rv=1, k=1, CV0=0.05, cal.sw=2)

}

\keyword{methods}
\keyword{survey} 