\name{quad_roots}
\alias{quad_roots}
\title{
Compute the roots of a quadratic equation   
}
\description{
Compute the roots of a quadratic equation
}
\usage{
quad_roots(a, b, c)
}
\arguments{
  \item{a}{coefficient of the quadratic term}
  \item{b}{coefficient of the linear term}
  \item{c}{coefficient of the constant term}
}

\details{
\code{quad_roots} computes the roots of a quadratic equation of the form \eqn{ax^2 + bx + c = 0}.
}
\value{vector with the two roots}

\examples{
    quad_roots(1, -8, 12)
}

\keyword{methods}
\keyword{survey} 