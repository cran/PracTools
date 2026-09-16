\name{pclass}
\alias{pclass}
\title{
Form nonresponse adjustment classes based on propensity scores
}
\description{
Fit a binary regression model for response probabilities and divide units into a specified number of classes.
}
\usage{
pclass(formula, data, link="logit", numcl=5, type, design=NULL, seed=NULL, rng.preds=0.05)
}
\arguments{
  \item{formula}{symbolic description of the binary regression model to be fitted as used in \code{glm}}
  \item{data}{an optional data frame containing columns for a response variable and covariates used in the model to predict response propensities; must be specified if \code{type="unwtd"}}
  \item{link}{a specification for the model link function; allowable values are \code{"logit"}, \code{"probit"}, or \code{"cloglog"}}
  \item{numcl}{number of classes into which units are split based on estimated propensities}
  \item{type}{whether an unweighted or weighted binary regression should be fit; allowable values are \code{"unwtd"} or \code{"wtd"}}
  \item{design}{sample design object; required if \code{type="wtd"}}
  \item{seed}{random number seed}
  \item{rng.preds}{cutoff value for the range of the predicted response probabilities}
}
\details{
A typical \code{formula} has the form \code{response ~ terms} where response is a two-level variable coded as 0 or 1, or is a factor where 
the first level denotes nonresponse and the second level is response. If \code{type="unwtd"}, 
\code{glm} is used to fit an unweighted regression. If \code{type="wtd"}, \code{svyglm} in the \code{survey} package is used to 
fit a survey-weighted regression. The \code{jitter} function is used to perturb regression predictions slightly. To reproduce the same set of weights in 
different runs, set \code{seed} to some non-null integer. Recommended seeds can be found in positions 3:626 of the built-in global vector 
\code{.Random.seed}.

If the range (defined as max - min) of the predicted response probabilities is less than \code{rng.preds}, the predictions are considered to have no 
meaningful variation and a single response adjustment class is returned. \code{rng.preds} is adjustable to accomodate different user preferences.
}
\value{
A list with components:
    \item{p.class}{propensity class for each unit}
    \item{propensities}{estimated response probability for each unit}
}

\references{
Valliant, R., Dever, J., Kreuter, F. (2018, chap. 13). \emph{Practical Tools for Designing and Weighting Survey Samples, 2nd edition}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}
\seealso{
\code{\link{NRadjClass}}
}
\examples{
    # classes based on unweighted logistic regression
require(PracTools)
data(nhis)
out <- pclass(formula = resp ~ age + as.factor(sex) + as.factor(hisp) + as.factor(race),
           data = nhis, type = "unwtd", link = "logit", numcl = 5, seed = 1387014860, 
           rng.preds = 0.05)
table(out$p.class, useNA="always")
summary(out$propensities)
    # classes based on survey-weighted logistic regression
require(survey)
nhis.dsgn <- svydesign(ids = ~psu, strata = ~stratum, data = nhis, nest = TRUE, weights = ~svywt)
out <- pclass(formula = resp ~ age + as.factor(sex) + as.factor(hisp) + as.factor(race),
           type = "wtd", design = nhis.dsgn, link="logit", numcl=5)
table(out$p.class, useNA="always")
summary(out$propensities)
    # degenerate case where all response propensities are the same
df <- data.frame(resp = rbinom(500, 1, 0.75), x = rnorm(500))
out1 <- pclass(resp ~ 1, data = df, type = "unwtd", numcl = 5, rng.preds = 0.05)    
}

\keyword{methods}
\keyword{survey} 
