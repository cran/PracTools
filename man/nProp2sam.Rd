\name{nProp2sam}
\alias{nProp2sam}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simple random sample size for difference in proportions
}
\description{
Compute a simple random sample size for estimating the difference in proportions when samples overlap
}
\usage{
nProp2sam(px, py, pxy, g, r, alt, sig.level=0.05, pow=0.80)
}
\arguments{
  \item{px}{
 proportion in group 1
}
  \item{py}{
 proportion in group 2
}
  \item{pxy}{
 proportion in the overlap has the characteristic in both samples
}
  \item{g}{
proportion of sample 1 that is in the overlap with sample 2
}
  \item{r}{
  ratio of the size of sample 1 to that of sample 2
}
  \item{alt}{
  should the test be 1-sided or 2-sided; allowable values are \code{alt="one.sided"} or  \code{alt="two.sided"}.
}
  \item{sig.level}{
  significance level of the hypothesis test
}
  \item{pow}{
  desired power of the test
}
}
\details{
\code{nProp2sam} computes sample sizes in two groups that are required for testing whether
the difference in group proportions is significant.
The power of the test is one of the input parameters.
The samples have a specified proportion of units in common.
}

\value{
 List with values:
  \item{n1}{sample size in group 1}
  \item{n2}{sample size in group 2}
  \item{px.py.pxy}{input values of the \code{px}, \code{py}, \code{pxy} parameters}
  \item{gamma}{proportion of sample 1 that is in the overlap with sample 2}
  \item{r}{ratio of the size of sample 1 to that of sample 2}
  \item{alt}{type of test: one-sided or two-sided}
  \item{sig.level}{significance level of test}
  \item{power}{power of the test}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2013, chap. 4). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.

Woodward, M. (1992). Formulas for Sample Size, Power, and Minimum Detectable Relative Risk in Medical Studies.  \emph{The Statistician}, 41, 185-196.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{nDep2sam}}
}
\examples{
nProp2sam(px=0.5, py=0.55, pxy=0.45, g=0.5, r=1, alt="two.sided")
}

\keyword{methods}
\keyword{survey} 