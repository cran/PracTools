\name{dub}
\alias{dub}

\title{
Sample sizes for a double sampling design
}
\description{
Compute samples sizes at each phase of a two-phase design where strata are created using the first phase.
}
\usage{
dub(c1, c2, Ctot, Nh, Sh, Yh.bar)
}

\arguments{
 \item{c1}{cost per unit in phase-1}
 \item{c2}{cost per unit in phase-2}
 \item{Ctot}{Total variable cost}
 \item{Nh}{Vector of stratum population counts or proportions}
 \item{Sh}{Vector of stratum population standard deviations}
 \item{Yh.bar}{Vector of stratum population means}
}

\details{
Compute the first and second phase sample sizes for a double sampling design. A first phase sample is selected by simple random sampling (\emph{srs}). Strata are formed based on information collected in the first phase. The Neyman allocation to strata of the second phase sample is computed ignoring costs. Optimal total sample sizes are computed for the first and second phases, given per-unit costs for the first and second phases and a fixed total budget for both phases combined.}

\value{
 A list object with elements:
  \item{V1}{Variance component associated with phase-1}
  \item{V2}{Variance component associated with phase-2}
  \item{n1}{Phase-1 sample size}
  \item{n2}{Total phase-2 sample across all strata}
  \item{"n2/n1"}{Fraction that phase-2 is of phase-1}
  \item{ney.alloc}{Vector of stratum sample sizes for phase-2 sample}
  \item{Vopt}{Variance of mean with the calculated phase-1 and phase-2 sample sizes}
  \item{nsrs}{Size of an \emph{srs} that has cost \code{Ctot}, assuming each unit costs \code{c2}}
  \item{Vsrs}{Variance of mean in an \emph{srs} of cost \code{Ctot}, assuming each unit costs \code{c2}}
  \item{Vratio}{Ratio of \code{Vopt} to \code{Vsrs}}
  \item{Ctot}{Input value of total cost}
  \item{cost.chk}{Computed value of phase-1 plus phase-2 sample with optimal sample sizes; should agree with Ctot}
}

\references{
Cochran W (1977, sect. 12.3) \emph{Sampling Techniques}. New York: John Wiley & Sons, Inc.

Neyman J (1938) Contribution to the theory of sampling human populations. \emph{Journal
of the American Statistical Association}, 33(201), 101-116.

Valliant, R., Dever, J., Kreuter, F. (2018, sect. 17.5.2).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\examples{
Wh <- rep(0.25,4)
Ph <- c(0.02,0.12,0.37,0.54)
Sh <- sqrt(Ph*(1-Ph))
c1 <- 10
c2 <- 50
Ctot <- 20000
dub(c1, c2, Ctot, Nh=Wh, Sh, Yh.bar=Ph)
}

\keyword{methods}
\keyword{survey} 