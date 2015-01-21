\name{NRadjClass}
\alias{NRadjClass}
\title{
Class-based nonresponse adjustments
}
\description{
Compute separate nonresponse adjustments in a set of classes.
}
\usage{
NRadjClass(ID, NRclass, resp, preds=NULL, wts=NULL, type)
}
\arguments{
  \item{ID}{identification value for a unit}
  \item{NRclass}{vector of classes to use for nonresponse adjustment. Length is number of respondents plus nonrespondents}
  \item{resp}{indicator for whether unit is anonrespondent (must be coded 0) or respondent (must be coded 1)}
  \item{preds}{response probabilities, typically estimated from a binary regression model as in \code{pclass}}
  \item{wts}{vector of survey weights, typically base weights or base weights adjusted for unknown eligibility}
  \item{type}{type of adjustment computed within each value of \code{NRclass}. Allowable codes are 1, 2, 3, 4, or 5.
      (1 = unweighted average of response propensities, i.e., \code{preds};
       2 = weighted average response propensity;
       3 = unweighted response rate;
       4 = weighted response rate;
       5 = median response propensity)
  }
}
\details{
The input vectors should include both respondents and nonrespondents in a sample. A single value between 0 and 1 is computed in each nonresponse adjustment class to be used as a nonresponse adjustment. Five alternatives are available for computing the adjustment based on the value of \code{type}. The value of the adjustment is merged with individual unit data and stored in the \code{RR} field of the output data frame.
}
\value{
A data frame of respondents only with four columns:
    \item{NRcl.no}{number of the nonresponse adjustment class for each unit}
    \item{ID}{identification value for a unit}
    \item{resp}{value of the \code{resp} variable (always 1)}
    \item{RR}{nonresponse adjustment for each unit}
}
\references{
Valliant, R., Dever, J., Kreuter, F. (2013, chap. 13).  \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\seealso{
\code{\link{pclass}}
}
\examples{
require(PracTools)
data(nhis)
out <- pclass(formula = resp ~ age + as.factor(sex) + as.factor(hisp) + as.factor(race),
         data = nhis, type = "unwtd", link="logit", numcl=5)
    # unweighted average of response propensities within each class
zz <- NRadjClass(ID=nhis[,"ID"], NRclass = as.numeric(out$p.class), resp=nhis[,"resp"],
            preds=out$propensities, wts=NULL, type=1)
}

\keyword{methods}
\keyword{survey} 