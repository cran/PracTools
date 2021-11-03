\name{strAlloc}
\alias{strAlloc}

\title{
Allocate a sample to strata
}
\description{
Compute the proportional, Neyman, cost-constrained, and variance-constrained
allocations in a stratified simple random sample.
}
\usage{
strAlloc(n.tot = NULL, Nh = NULL, Sh = NULL, cost = NULL, ch = NULL,
     V0 = NULL, CV0 = NULL, ybarU = NULL, alloc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.tot}{
fixed total sample size
}
  \item{Nh}{
vector of population stratum sizes (\eqn{N_h}) or pop stratum proportions (\eqn{W_h})
}
  \item{Sh}{
stratum unit standard deviations (\eqn{S_h}), required unless \code{alloc = "prop"}
}
  \item{cost}{
total variable cost
}
  \item{ch}{
vector of costs per unit in stratum \emph{h} \eqn{(c_h)}
}
  \item{V0}{
fixed variance target for estimated mean
}
  \item{CV0}{
fixed CV target for estimated mean
}
  \item{ybarU}{
population mean of \emph{y} (\eqn{\bar{y}_U})
}
  \item{alloc}{
type of allocation; must be one of \code{"prop"}, \code{"neyman"}, \code{"totcost"}, \code{"totvar"}
}
}
\details{
\code{alloc="prop"} computes the proportional allocation of a fixed total sample size, \code{n.tot}, to the strata. \code{alloc="neyman"} computes the allocation of a fixed total sample size, \code{n.tot}, to the strata that minimizes the variance of an estimated mean. \code{alloc="totcost"} computes the allocation of a fixed total sample size, \code{n.tot}, to the strata that minimizes the variance of an estimated mean subject to the fixed total \code{cost}.  \code{alloc="totvar"} computes the allocation that minimizes total cost subject to the target coefficient of variation, \code{CV0}, or the target variance, \code{V0}, of the estimated mean.
}

\value{
    For proportional allocation, a list with values:
    \item{alloc}{type of allocation: \code{"prop"}, \code{"neyman"}, \code{"totcost"}, \code{"totvar"}}
    \item{Nh}{vector of population sizes (\eqn{N_h}) or pop stratum proportions (\eqn{W_h})}
    \item{nh}{vector of stratum sample sizes}
    \item{"nh/n"}{proportion of sample allocated to each stratum}
    For other allocations, the three components above plus:
    \item{Sh}{stratum unit standard deviations (\eqn{S_h})}
    \item{"anticipated SE of estimated mean"}{Anticipated SE of the estimated mean for the computed allocation}
}
\references{
 Cochran, W.G. (1977). \emph{Sampling Techniques}. John Wiley & Sons.

 Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. Springer.
}
\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\seealso{
\code{\link{nCont}}, \code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
}

\examples{
# Neyman allocation
Nh <- c(215, 65, 252, 50, 149, 144)
Sh <- c(26787207, 10645109, 6909676, 11085034, 9817762, 44553355)
strAlloc(n.tot = 100, Nh = Nh, Sh = Sh, alloc = "neyman")

# cost constrained allocation
ch <- c(1400, 200, 300, 600, 450, 1000)
strAlloc(Nh = Nh, Sh = Sh, cost = 100000, ch = ch, alloc = "totcost")

# allocation with CV target of 0.05
strAlloc(Nh = Nh, Sh = Sh, CV0 = 0.05, ch = ch, ybarU = 11664181, alloc = "totvar")
}

\keyword{methods}
\keyword{survey}
