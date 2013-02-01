\name{HMT}
\alias{HMT}

\title{
Generate an HMT population
}
\description{
Generate a population that follows the model in Hansen, Madow, and Tepping (1983)
}
\usage{
HMT(N=5000, H=10)
}

\arguments{
 \item{N}{population size}
 \item{H}{number of strata}
}

\details{
\code{HMT} generates a population based on the model: \eqn{E(y)= \alpha + \beta x}, \eqn{Var(y)=\sigma^2x^{3/2}}. Both \emph{x} and \emph{y} have gamma distributions. Strata are formed to have approximately the same total of \eqn{x}.}

\value{
 N x 3 matrix with columns:
  \item{strat}{stratum ID}
  \item{x}{auxiliary variable \emph{x}}
  \item{y}{analysis variable \emph{y}}
}

\references{
Hansen,M.H., Madow,W.G., and Tepping,B.J. (1983). An evaluation of model-dependent and probability sampling inferences in sample surveys. \emph{Journal of the American Statistical Association}, 78, 776-793.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\examples{
    # generate HMT population with 1000 units and 5 strata and plot results
pop <- HMT(N=1000, H=5)
plot(pop[,"x"],pop[,"y"])
}

\keyword{methods}
\keyword{survey} 