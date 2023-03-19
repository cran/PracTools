\name{nCont}
\alias{nCont}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a simple random sample size for an estimated mean
}
\description{
Compute a simple random sample size using either a target coefficient of variation, \eqn{CV_0}, or target variance, \eqn{V_0}, for an estimated mean.
}
\usage{
nCont(CV0=NULL, V0=NULL, S2=NULL, ybarU=NULL, N=Inf, CVpop=NULL)
}
\arguments{
  \item{CV0}{
 target value of coefficient of variation of \eqn{\bar{y_s}}
}
  \item{V0}{
target value of variance of \eqn{\bar{y_s}}
}
  \item{S2}{
 unit (population) variance
}
  \item{ybarU}{
 population mean of target variable
}
  \item{N}{
number of units in finite population
}
  \item{CVpop}{
  unit (population) coefficient of variation
}
}
\details{
If \eqn{CV_0} is the desired target, then the unit CV,  \code{CVpop},
or the population mean and variance, \code{ybarU} and \code{S2}, must also be provided.
If \eqn{V_0} is the constrained value, then \code{S2} must be also be included in the function call.
}
\value{
  numeric sample size
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{nLogOdds}}, \code{\link{nProp}}, \code{\link{nPropMoe}}, \code{\link{nWilson}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nCont(CV0=0.05, CVpop=2)
nCont(CV0=0.05, CVpop=2, N=500)
nCont(CV0=0.10/1.645, CVpop=1)

    # Compute sample size for a ratio estimator in smho98 population
    # Identify large units to select with certainty first
data(smho98)
cert <- smho98[,"BEDS"] > 2000
tmp <- smho98[!cert, ]
tmp <- tmp[tmp[, "BEDS"] > 0, ]

x <- tmp[,"BEDS"]
y <- tmp[, "EXPTOTAL"]
m <- lm(y ~ 0 + x, weights = 1/x)
ybarU <- mean(y)
S2R <- sum(m$residuals^2/(length(x)-1))
nCont(CV0=0.15, S2=S2R, ybarU=ybarU)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{methods}
\keyword{survey} 