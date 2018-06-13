\name{clusOpt3fixedPSU}
\alias{clusOpt3fixedPSU}
\title{
Compute optimal number of sample secondary sampling units (SSUs) and elements per SSU for a fixed set of primary sampling units (PSUs) in a three-stage sample
}
\description{
Compute the sample sizes that minimize the variance of the \emph{pwr}-estimator of a total in a three-stage sample when the PSU sample is fixed.
}
\usage{
clusOpt3fixedPSU(unit.cost, m, delta1, delta2, unit.rv, k1=1, k2=1, CV0=NULL,
         tot.cost=NULL, cal.sw)
}
\arguments{
   \item{unit.cost}{3-vector of unit costs:
      \code{C1} = unit cost per PSU;
      \code{C2} = unit cost per SSU;
      \code{C3} = unit cost per element
   }
   \item{m}{number of sample PSUs (fixed)}
   \item{delta1}{homogeneity measure among elements within PSUs}
   \item{delta2}{homogeneity measure among elements within SSUs}
   \item{unit.rv}{unit relvariance}
   \item{k1}{ratio of \eqn{B^2 + W^2} to unit relvariance}
   \item{k2}{ratio of \eqn{W_{2}^2 + W_{3}^2} to unit relvariance}
   \item{CV0}{target CV}
   \item{tot.cost}{total budget for variable costs, including PSU costs}
   \item{cal.sw}{specify type of optimum:
               1 = find optimal \code{m.opt} for fixed total budget;
               2 = find optimal \code{m.opt} for target CV0}
}

\details{
\code{clusOpt3} will compute \eqn{\bar{n}_{opt}} and \eqn{\bar{\bar{q}}_{opt}} for a three-stage
sample which uses simple random sampling at each stage or \emph{ppswr} at the first stage and \emph{srs} at the second and third stages. The set of sample PSUs is assumed to be fixed.
"Variable costs" in \code{tot.cost} includes the budget for all costs that vary with the number
of sample PSUs, SSUs, and elements, i.e., \eqn{C_{1}m + C_{2}m\bar{n} + C_{3}m\bar{n}\bar{\bar{q}}}.
}
\value{
 List with values:
  \item{C1}{unit cost per PSU}
  \item{C2}{unit cost per SSU}
  \item{C3}{unit cost per element}
  \item{m}{ number of sample PSUs (fixed)}
  \item{delta1}{homogeneity measure among elements within PSUs}
  \item{delta2}{homogeneity measure among elements within SSUs}
  \item{unit relvar}{unit relvariance}
  \item{k1}{ratio of \eqn{B^2 + W^2} to unit relvariance}
  \item{k2}{ratio of \eqn{W_{2}^2 + W_{3}^2} to unit relvariance}
  \item{cost}{budget constraint, \code{tot.cost} if \code{cal.sw}=1; computed cost if \code{cal.sw}=2}
  \item{n}{optimum number of sample SSUs per PSU}
  \item{q}{optimum number of sample elements per SSU}
  \item{CV}{achieved CV, used if \code{cal.sw}=1; or target CV, used if \code{cal.sw}=2}
  \item{CV check}{computed CV based on optimal sample sizes; used only if \code{cal.sw}=2}
}

\references{
Hansen,M.H., Hurwitz,W.N., and Madow,W.G. (1953, p. 225).  \emph{Sample Survey Methods and Theory}, Vol.II. John Wiley & Sons.

Valliant, R., Dever, J., Kreuter, F. (2013, sect. 9.3.2).  \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{clusOpt2}}, \code{\link{clusOpt2fixedPSU}}, \code{\link{clusOpt3}}
}

\examples{
    # optima for a fixed total budget
clusOpt3fixedPSU(unit.cost=c(500, 100, 120), m=100, delta1=0.01, delta2=0.05, unit.rv=1,
     k1=1, k2=1, tot.cost=500000,cal.sw=1)
    # optima for a target CV
clusOpt3fixedPSU(unit.cost=c(500, 100, 120), m=100, delta1=0.01, delta2=0.05, unit.rv=1,
     k1=1, k2=1, CV0=0.05,cal.sw=2)
}

\keyword{methods}
\keyword{survey} 