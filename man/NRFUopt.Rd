\name{NRFUopt}
\alias{NRFUopt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sample sizes for a nonresponse follow-up study
}
\description{
Compute optimal values of the first-phase sample size and the second-phase sampling fraction in a two-phase sample.
}
\usage{
NRFUopt(Ctot=NULL, c1, c2, theta, CV0=NULL, CVpop=NULL, N=Inf, type.sw)
}

\arguments{
 \item{Ctot}{total variable cost}
 \item{c1}{cost per unit in phase-1}
 \item{c2}{cost per unit in phase-2}
 \item{theta}{probability of response for each unit}
 \item{CV0}{target coefficient of variation for the estimated total or mean}
 \item{CVpop}{Unit coefficient of variation}
 \item{N}{Population size; default is \code{Inf}}
 \item{type.sw}{type of allocation; \code{"cost"} = target total variable cost, \code{"cv"} = target coefficient of variation}
}

\details{
\code{NRFUopt} computes the optimal values of the first-phase sample size and the second-phase sampling fraction in a two-phase sample. Both stages are assumed to be selected using simple random sampling without replacement. If \code{type.sw="cost"}, the optima are computed for a target total, expected cost across both phases.If \code{type.sw="cv"}, the optima are computed for a target coefficient of variation for an estimated mean.
}

\value{
 List object with values:
  \item{allocation}{type of allocation: either "fixed cost" or "fixed CV"}
  \item{"Total variable cost"}{expected total cost: fixed cost if \code{type.sw="cost"} or computed cost if \code{type.sw="cv"}}
  \item{"Response rate"}{first-phase response rate}
  \item{CV}{anticipated coefficient of variation (CV) if \code{type.sw="cost"} or target CV if \code{type.sw="cv"}}
  \item{v.opt}{optimal fraction of first-phase nonrespondents to select for second-phase follow-up}
  \item{n1.opt}{optimal number of units to sample at first-phase}
  \item{"Expected n2"}{expected number of respondents obtained at second-phase}
  \item{"Expected total cases (2-phase)"}{expected number of respondents across both phases}
  \item{"srs sample for same cv"}{size of single-phase simple random sample (\emph{srs}) needed to obtain same CV as the two-phase sample}
  \item{"Cost Ratio: Two phase to srs"}{ratio of expected cost for two-phase sample to cost of single-phase \emph{srs}}
}

\references{
Saerndal,C.E., Swensson, B., and Wretman, J. (1992, examples 15.4.4 and 15.4.5). \emph{Model Assisted Survey Sampling}. New York: Springer.

Valliant, R., Dever, J., Kreuter, F. (2013, chap.17). \emph{Practical Tools for Designing and Weighting Survey Samples}. New York: Springer.
}

\author{
Richard Valliant, Jill A. Dever, Frauke Kreuter
}

\examples{
    # optima for fixed target CV
NRFUopt(Ctot=NULL, c1=50, c2=200, theta=0.5, CV0=0.05, CVpop=1, type.sw = "cv")
    # optima for fixed total cost
NRFUopt(Ctot=100000, c1=50, c2=200, theta=0.5, CV0=NULL, CVpop=1, type.sw = "cost")

}

\keyword{methods}
\keyword{survey} 