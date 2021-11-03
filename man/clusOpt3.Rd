\name{clusOpt3}
\alias{clusOpt3}
\title{
Compute optimal sample sizes for a three-stage sample
}
\description{
Compute the sample sizes that minimize the variance of the \emph{pwr}-estimator of a total in a three-stage sample.
}
\usage{
clusOpt3(unit.cost, delta1, delta2, unit.rv, k1=1, k2=1, CV0=NULL, tot.cost=NULL, cal.sw)
}
\arguments{
   \item{unit.cost}{vector with three components for unit costs:
	   \code{C1} = unit cost per primary sampling unit (PSU);
	   \code{C2} = unit cost per  secondary sampling unit (SSU);
	   \code{C3} = unit cost per element
    }
   \item{delta1}{homogeneity measure among elements within PSUs}
   \item{delta2}{homogeneity measure among elements within SSUs}
   \item{unit.rv}{population unit relvariance}
   \item{k1}{ratio of \eqn{B^2 + W^2} to the population unit relvariance}
   \item{k2}{ratio of \eqn{W_{2}^2 + W_{3}^2} to the population unit relvariance}
   \item{CV0}{target CV}
   \item{tot.cost}{total budget for variable costs}
   \item{cal.sw}{specify type of optimum:
               1 = find optimal \code{m.opt} for fixed total budget;
               2 = find optimal \code{m.opt} for target CV0}
}

\details{
\code{clusOpt3} will compute \eqn{m_{opt}}, \eqn{\bar{n}_{opt}}, and \eqn{\bar{\bar{q}}_{opt}}
for a three-stage sample which uses simple random sampling at each stage or \emph{ppswr}
at the first stage and \emph{srs} at the second and third stages.
}
\value{
 List with values:
  \item{C1}{unit cost per PSU}
  \item{C2}{unit cost per SSU}
  \item{C3}{unit cost per element}
  \item{delta1}{homogeneity measure among elements within PSUs}
  \item{delta2}{homogeneity measure among elements within SSUs}
  \item{unit relvar}{unit relvariance}
  \item{k1}{ratio of \eqn{B^2 + W^2} to the population unit relvariance}
  \item{k2}{ratio of \eqn{W_{2}^2 + W_{3}^2} to the population unit relvariance}
  \item{cost}{total budget for variable costs if \code{cal.sw}=1; or computed cost if \code{cal.sw}=2}
  \item{m.opt}{optimum number of sample PSUs}
  \item{n.opt}{optimum number of sample SSUs per PSU}
  \item{q.opt}{optimum number of sample elements per SSU}
  \item{CV}{achieved CV if \code{cal.sw}=1 or target CV if \code{cal.sw}=2}
}

\references{
Hansen,M.H., Hurwitz,W.N., and Madow,W.G. (1953, p. 225).  \emph{Sample Survey Methods and Theory}, Vol.II. John Wiley & Sons.

Valliant, R., Dever, J., Kreuter, F. (2018, sect. 9.3.2).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{clusOpt2}}, \code{\link{clusOpt2fixedPSU}}, \code{\link{clusOpt3fixedPSU}}
}

\examples{
    # optima for a fixed total budget
clusOpt3(unit.cost=c(500, 100, 120), delta1=0.01, delta2=0.10, unit.rv=1,
       k1=1, k2=1, tot.cost=100000,cal.sw=1)
    # optima for a target CV
clusOpt3(unit.cost=c(500, 100, 120), delta1=0.01, delta2=0.10, unit.rv=1,
       k1=1, k2=1, CV0=0.01,cal.sw=2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{methods}
\keyword{survey} 