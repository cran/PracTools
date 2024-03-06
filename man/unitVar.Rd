\name{unitVar}
\alias{unitVar}
\title{
Compute the unit (population) variance for a variable
}
\description{
Compute the unit (population) variance for a variable based on either a full population file or a sample from a finite population.
}
\usage{
unitVar(pop.sw = NULL, w = NULL, p = NULL, y = NULL)
}
\arguments{
  \item{pop.sw}{
    TRUE if the full population is input; FALSE if a sample is input
}
  \item{w}{
    vector of sample weights if \code{y} is a sample; used only if \code{pop.sw = FALSE}
}
  \item{p}{
    vector of 1-draw selection probabilities; optionally provided if \code{pop.sw = TRUE}
}
  \item{y}{
     vector of values of an analysis variable; must be numeric
}
}
\details{
    \code{unitVar} computes unit (population) variances of an analysis variable \eqn{y} from either a population or a sample. \code{S2} is the unweighted population variance, \eqn{S^2 = \sum_{i \in U}(y_i - \bar{y}_U)^2/(N-1)} where \eqn{U} is the universe of elements, \eqn{N} is the population size, and \eqn{\bar{y}_U} is the population mean. If the input is a sample, \code{S2} is estimated as \eqn{\hat{S}^2 = (n/(n-1))\sum_{i \in s} w_i(y_i - \bar{y}_w)^2/(\sum_{i \in s} w_i)} where \eqn{s} is the set of sample elements, \eqn{n} is the sample size, and \eqn{\bar{y}_w} is the weighted sample mean.
    
    \code{V1} is a weighted population variance used in calculations for samples where elements are selected with varying probabilities. If the \eqn{y} is a population vector, \eqn{V_1 = \sum_U p_i(y_i/p_i - t_U)^2} where \eqn{p_i} is the 1-draw probability for element \eqn{i} and \eqn{t_U} is the population total of \eqn{y}. If \eqn{y} is for a sample, \eqn{\hat{V}_1 = \sum_s (y_i/p_i - n^{-1}\sum_k y_k/p_k)^2 / (n-1)} with \eqn{p_i} computed as \eqn{1/(n w_i)}.
}
\value{
A list with three or four components:
    \item{Note}{Describes whether output was computed from a full population or estimated from a sample.}
    \item{Pop size N}{Size of the population; included if \code{y} is for the full population.}
    \item{S2}{Unit variance of \code{y}; if \code{pop.sw = TRUE}, \code{S2} is computed from the full population; if \code{pop.sw = FALSE}, \code{S2} is estimated from the sample using the \code{w} weights.}
    \item{V1}{Population variance of \code{y} appropriate for a sample selected with varying probabilities; see Valliant, Dever, and Kreuter (VDK; 2018, sec. 3.4). If \code{pop.sw = TRUE} and \code{p} is provided, \code{V1} is computed with equation (3.32) in VDK. If \code{pop.sw = FALSE}, \code{V1} is estimated with equation (3.41) in VDK.}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 3). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}
\seealso{
    \code{\link{nCont}}, \code{\link{nContMoe}}, \code{\link{nContOpt}}, \code{\link{nPPS}}
}
\author{
Richard Valliant
}

\examples{
library(PracTools)
data("smho.N874")
y <- smho.N874[,"EXPTOTAL"]
x <- smho.N874[, "BEDS"]
y <- y[x>0]
x <- x[x>0]
pik <- x/sum(x)
require(sampling)
n <- 50
sam <- UPrandomsystematic(n * pik)
wts <- 1/(n*pik[sam==1])
unitVar(pop.sw = TRUE, w = NULL, p = pik, y=y)
unitVar(pop.sw = FALSE, w = wts, p = NULL, y=y[sam==1])
}

\keyword{methods}
\keyword{survey} 