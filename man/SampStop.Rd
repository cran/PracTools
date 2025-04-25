\name{SampStop}
\alias{SampStop}

\title{
Stopping rule for surveys
}
\description{
Compute the probability that continuing data collection will lead to a change in the value of an estimated mean.
}
\usage{
SampStop(lm.obj, formula, n1.data, yvar, n2.data, p = NULL, delta = NULL, seed = NULL)
}

\arguments{
  \item{lm.obj}{
  object of class \code{lm} from a regression predicting \eqn{y} based on \code{n1.data}
}
  \item{formula}{
  righthand side of the formula in \code{lm.obj}; it excludes the dependent variable \eqn{y}; no quotes are used.
}
  \item{n1.data}{
  data frame containing units in the part of the sample that has been completed; includes \eqn{y} and the covariates in \code{formula}.
}
  \item{yvar}{
  name or number of column in \code{n1.data} containing \eqn{y}.
}
  \item{n2.data}{
  data frame containing units in the part of the sample that is yet to be completed; includes only covariates in \code{formula}.
}
  \item{p}{
  Vector of anticipated response probabilities for the n2 sample; 0 < \code{p} < 1.
}
  \item{delta}{
  vector of potential differences in the estimated means for the n1 and n2 samples.
}
  \item{seed}{
  random number seed for selecting sample from incomplete cases.
}

}

\details{
\code{SampStop} allows an evaluation to be made of whether data collection can be stopped, without substantially affecting the value
  of an estimated mean, prior to completing collection for all units. Suppose that a sample of size \eqn{n} is divided between 
  the \eqn{n_1} units whose collection has been completed and the remaining \eqn{n_2 = n - n_1} units that are yet to be completed. The 
  function computes \eqn{Pr(|e_1 - e_2| < \delta)} where \eqn{e_1 - e_2} is the potential difference (\code{delta}) between the estimated mean 
  based on the completed sample and the estimated mean for the full sample if all units were to be completed. For \eqn{e_1} the mean is estimated 
  after imputing the \eqn{y}'s for the \eqn{n_2} incomplete units. The estimated mean \eqn{e_2} is computed assuming that an additional 
  \eqn{n_2 * p} units are completed, and the \eqn{y}'s for the remaining \eqn{n_2 - n_2*p} incomplete units are imputed. Estimating the variance 
  of \eqn{e_1 - e_2}  involves selecting a sample from \code{n2.data} using the random number seed in \code{seed}.
  
  The parameter \code{p} is the response rate that is anticipated for the \eqn{n_2} uncompleted units. The usual situation is that 
  there is some uncertainty about \code{p} which can be accounted for by inputting a vector of \code{p}'s. \eqn{\delta} is a difference in estimates 
  that, if not exceeded, would lead to stopping data collection. For an acceptably small value of \code{delta}, if 
  \eqn{Pr(|e_1 - e_2| < \delta)} is large enough, the decision can be made to stop data collection. The variable \eqn{y} in \code{yvar} is assumed to follow the linear model
  in \code{lm.obj}. A model with independent errors (or a simple random sample) is assumed for calculations.
}

\value{Matrix with \code{length(p)*length{delta}} rows and columns:
  \item{Pr(response)}{Probability of response by each of the remaining \eqn{n_2} cases}
  \item{Exp no. resps}{Expected number of respondents among the remaining \eqn{n_2} cases}, i.e. \eqn{n_2*p}
  \item{y1 mean}{Mean of the \eqn{n_1} respondents}
  \item{diff in means}{Value of the input parameter \code{delta}}
  \item{se of diff}{Standard error of the difference \code{delta}}
  \item{z-score}{Z-score for computing \eqn{Pr(|e_1 - e_2| < \delta})}
  \item{Pr(smaller diff)}{\eqn{Pr(|e_1 - e_2| < \delta)} for the inputs of \code{p} and \code{delta}}
}

\references{
Wagner, J. and Raghunathan, T. (2010). A new stopping rule for surveys. \emph{Statistics in Medicine}, 29(9), 1014-1024.
}

\author{
George Zipf, Richard Valliant
}

\examples{
library(PracTools)
    # Model with quantitative covariates
data(hospital)
HOSP <- hospital
HOSP$sqrt.x <- sqrt(HOSP$x)
sam   <- sample(nrow(HOSP), 50)
N1       <- HOSP[sam, ]
N2       <- HOSP[-sam, ]
    ## Create lm object using "known" data; no intercept model
lm.obj  <- lm(y ~ 0 + sqrt.x + x, data = N1)
del <- mean(HOSP$y) - mean(HOSP$y) * seq(.6, 1, by=0.05)
SampStop(lm.obj  = lm.obj,
                    formula = ~ 0 + sqrt.x + x,
                    n1.data = N1,
                    yvar    = "y",
                    n2.data = N2,
                    p       = seq(0.2, 0.6, by=0.05),
                    delta   = del,
                    seed = .Random.seed[413]) 
    # Model with factors
data(labor)
sam   <- sample(nrow(labor), 50)
n1.vars <- c("WklyWage", "HoursPerWk", "agecat", "sex")
n2.vars <- c("HoursPerWk", "agecat", "sex")
N1       <- labor[sam, n1.vars]
N2       <- labor[-sam, n2.vars]
lm.obj  <- lm(WklyWage ~ HoursPerWk + as.factor(agecat) + as.factor(sex), data = labor)
del <- mean(N1$WklyWage) - mean(N1$WklyWage) * seq(.75, .95, by=0.05)
result <- SampStop(lm.obj  = lm.obj,
                    formula = ~ HoursPerWk + as.factor(agecat) + as.factor(sex),
                    n1.data = N1,
                    yvar    = "WklyWage",
                    n2.data = N2,
                    p       = seq(0.2, 0.4, by=0.05),
                    delta   = del,
                    seed = .Random.seed[78]) 

p.nresp <- paste(result[,1], result[,2], sep=", ")
library(ggplot2)
ggplot2::ggplot(result, aes(result[,4], result[,7], colour = factor(p.nresp) )) +
  geom_point() +
  geom_line(linewidth=1.1) +
  labs(x = "delta", y = "Pr(|e1-e2|<= delta)", colour = "Pr(resp), n.resp")                      
}

\keyword{methods}
\keyword{survey} 
