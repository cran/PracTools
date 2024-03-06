\name{nDep2sam}
\alias{nDep2sam}
\title{
Simple random sample size for difference in means
}
\description{
Compute a simple random sample size for estimating the difference in means when samples overlap
}
\usage{
nDep2sam(S2x, S2y, g, r, rho, alt, del, sig.level=0.05, pow=0.80)
}
\arguments{
  \item{S2x}{
 unit variance of analysis variable \emph{x} in sample 1
}
  \item{S2y}{
 unit variance of analysis variable \emph{y} in sample 2
}
  \item{g}{
proportion of sample 1 that is in the overlap with sample 2
}
  \item{r}{
  ratio of the size of sample 1 to that of sample 2
}
  \item{rho}{
  unit-level correlation between \emph{x} and \emph{y}
}
  \item{alt}{
  should the test be 1-sided or 2-sided; allowable values are \code{alt="one.sided"} or  \code{alt="two.sided"}.
}
  \item{del}{
  size of the difference between the means to be detected
}
  \item{sig.level}{
  significance level of the hypothesis test
}
  \item{pow}{
  desired power of the test
}
}
\details{
\code{nDep2sam} computes sample sizes in two groups that are required for testing whether the difference in group means is significant.  The power of the test is one of the input parameters. The samples have a specified proportion of units in common. Both samples are assumed to be selected via simple random sampling.
}

\value{
 List with values:
  \item{n1}{sample size in group 1}
  \item{n2}{sample size in group 2}
  \item{S2x.S2y}{unit variances in groups 1 and 2}
  \item{delta}{difference in group means to be detected}
  \item{gamma}{proportion of sample 1 that is in the overlap with sample 2}
  \item{r}{ratio of the size of sample 1 to that of sample 2}
  \item{rho}{unit-level correlation between analysis variables in groups 1 and 2}
  \item{alt}{type of test: one-sided or two-sided}
  \item{sig.level}{significance level of test}
  \item{power}{power of the test}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 4).  \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.

Woodward, M. (1992). Formulas for Sample Size, Power, and Minimum Detectable Relative Risk in Medical Studies.  \emph{The Statistician}, 41, 185-196.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{nProp2sam}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
nDep2sam(S2x=200, S2y=200,
            g=0.75, r=1, rho=0.9,
            alt="one.sided", del=5,
            sig.level=0.05, pow=0.80)
}

\keyword{methods}
\keyword{survey} 