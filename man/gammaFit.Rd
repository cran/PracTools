\name{gammaFit}
\alias{gammaFit}

\title{
Iteratively estimate variance model parameter \eqn{\gamma}
}
\description{
Iteratively computes estimate of \eqn{\gamma} in a model with \eqn{E_M(y)=x^T\beta} and
\eqn{Var_M(y)=\sigma^2x^\gamma}.
}
\usage{
gammaFit(X, x, y, maxiter = 100, show.iter = FALSE, tol = 0.001)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
matrix of predictors in the linear model for \emph{y}
}
  \item{x}{
vector of \emph{x}'s for individual units in the assumed specification of \eqn{Var_M(y)}
}
  \item{y}{
vector of dependent variables for individual units
}
  \item{maxiter}{
maximum number of iterations allowed
}
  \item{show.iter}{
should values of \eqn{\gamma} be printed of each iteration? \code{TRUE} or \code{FALSE}
}
  \item{tol}{
size of relative difference in \eqn{\hat{\gamma}}'s between consecutive iterations
used to determine convergence. Algorithm terminates when relative difference
is less than \code{tol}.
}
}
\details{
The function \code{gammaFit} estimates the power \eqn{\gamma} in a model where the variance
of the errors is proportional to \eqn{x^\gamma} for some covariate \emph{x}.
Values of \eqn{\gamma} are typically in [0,2]. The function calls \code{\link{gamEst}}.
}
\value{
A list with the components:
\item{g.hat}{estimate of \eqn{\gamma} when iterative procedure stopped}
\item{converged}{\code{TRUE} or \code{FALSE} depending on whether convergence was obtained}
\item{steps}{number of steps used by the algorithm}
}
\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}
\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{gamEst}}
}
\examples{
data(hospital)
x <- hospital$x
y <- hospital$y

X <- cbind(sqrt(x), x)
gammaFit(X = X, x = x, y = y, maxiter=100, tol=0.002)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

