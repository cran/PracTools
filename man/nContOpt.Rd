\name{nContOpt}
\alias{nContOpt}                           
\title{      
Compute the sample size required to estimate the mean of a continuous variable by optimizing the numbers of take-alls and non-take-all units selected by probability sampling
}
\description{
Compute a sample size required to achieve a precision target for an estimated mean, \eqn{\hat{\bar{y}}_s}, based on splitting the sample between take-alls and non-take-alls. The sample design for non-take-alls can be either simple random sampling or probability proportional to size sampling.
}
\usage{
nContOpt(X, Y = NULL, CV0 = NULL, V0 = NULL, design = NULL)
}
\arguments{
  \item{X}{
      population variable used for determining take-all cutoff and for selecting a probability proportional to size sample if \code{design = "PPS"}. \code{X} is a vector that contains a value for every unit in the population.
}
  \item{Y}{
       variable used for computing a population variance; required if \code{design = "PPS"}. \code{Y} is ignored if  \code{design = "SRS"}. \code{X} is a vector that contains a value for every unit in the population and is the same length as \code{Y}.
}  
  \item{CV0}{
      target value of coefficient of variation of \eqn{\hat{\bar{y}}_s}
}
  \item{V0}{
     target value of variance of \eqn{\hat{\bar{y}}_s}
}
  \item{design}{
     Sample design to be used for non-take-alls; must be either \code{"SRS"} or \code{"PPS"}.
}
}
\details{
Compute a sample size based on splitting the sample between take-alls and non-take-alls in a way that achieves either a target coefficient of variation or a target variance for an estimated mean. The function sorts the file in descending order by \code{X} and then systematically designates units as  take-alls (certainty selections) starting from largest to smallest, and computes the sample size of non-take-alls needed to achieve the precision target. Initially, no unit in the ordered list is a certainty, and if \code{design = "SRS"}, the first value in \code{nContOpt.curve} is the same as \code{nCont} produces under identical inputs. In each pass, the algorithm increases the number of certainties. In the second pass, the first value is taken as a certainty and the non-take-all sample size is based on units 2:N, where N is the population size. On the third pass, the first two values are taken as certainties and the non-take-all sample size is based on units 3:N. The function cycles through units 1:(N-1) with take-alls increasing by 1 each cycle, and determines the minimum total sample size needed to achieve the specified precision target. The optimum sample size \code{nContOpt.n} combines certainties and non-certainties for its value.

The sample design can be either simple random sampling or probability proportional to size sampling. When \code{design = "SRS"}, calculations are based only on \code{X}. The SRS variance formula is for without replacement sampling so that a finite population correction factor (\emph{fpc}) is included. When \code{design = "PPS"}, \code{X} is used for the measure of size and \code{Y} is the variable for computing the variance used to determine the sample size. The PPS variance is computed for a with-replacement design, but an ad hoc \emph{fpc} is included. Either \code{CV0} or \code{V0} must be provided but not both.
}
\value{
 A list with five components:
   \item{nContOpt.Curve}{The sample size for the given inputs based on the number of take-alls incrementing from 1 to N-1}
   \item{Take.alls}{A TRUE/FALSE vector for whether the element in the \code{X} vector is a take-all}
   \item{nContOpt.n}{The minimum sample size (take-alls + non-take-alls) required for the given inputs, rounded to 4 decimal places}
   \item{Min.Takeall.Val}{The minimum value of \code{X} for the take-alls}
   \item{n.Take.all}{The number of take-all units in the optimal sample}
}
\author{George Zipf, Richard Valliant}
\seealso{
\code{\link{nCont}}, \code{\link{nContMoe}}
}
\examples{
nContOpt(X = TPV$Total.Pot.Value, CV0 = 0.05, design = "SRS")
nContOpt(X = TPV$Total.Pot.Value, V0  = 5e+14, design = "SRS")
g <- nContOpt(X = TPV$Total.Pot.Value, CV0 = 0.05, design = "SRS")
plot(g$nContOpt.Curve,
     type = "o",
     main = "Sample Size Curve",
     xlab = "Take-all / Sample Split Starting Value",
     ylab = "Total sample size (take-alls + non-tale-alls)" )
nContOpt(X = TPV$Total.Pot.Value, Y = TPV$Y, CV0 = 0.05, design = "PPS")
}
\keyword{methods}
\keyword{survey} 