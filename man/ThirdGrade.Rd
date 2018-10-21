\name{ThirdGrade}
\alias{ThirdGrade}
\docType{data}
\title{
Third grade population
}
\description{
The \code{ThirdGrade} data file is a population of students who participated in the Third International Mathematics and Science Study (TIMSS).
}

\usage{data(ThirdGrade)}

\format{
  A data frame with 2,427 students on the following variables:
  \describe{
  \item{\code{region}}{Geographic region of the U.S.
                  (1 = Northeast;
                  2 = South;
                  3 = Central;
                  4 = West)
  }
  \item{\code{school.id}}{School identifier (1 - 135)}
  \item{\code{student.id}}{Student identifier (1 - 2427)}
  \item{\code{sex}}{Sex of student
                  (1 = female;
                  2 = male)
  }
  \item{\code{language}}{Is language of test spoken at home?
                  (1 = always;
                  2 = sometimes;
                  3 = never)
  }
  \item{\code{math}}{Mathematics test score}
  \item{\code{ethnicity}}{Ethnicity of student
                  (1 = White, non-Hispanic;
                  2 = Black;
                  3 = Hispanic;
                  4 = Asian;
                  5 = Native American;
                  6 = Other)
  }
  \item{\code{science}}{Science test score}
  \item{\code{community}}{Type of location of school
                  (2 = village or rural area;
                  3 = outskirts of a town or city; 
                  4 = close to center of a town or city)
   }
   \item{\code{enrollment}}{Number of students in entire school}
 }
}
\details{
The Third Grade population consists of 2,427 students in the U.S. who participated in the Third International Mathematics and Science Study (Caslyn, Gonzales, Frase 1999).  The methods used in conducting the original study are given in TIMSS International Study Center (1996). Clusters are schools while units within clusters are the students.  }
\source{
    TIMSS International Study Center 1996.
}
\references{
    Caslyn, C., Gonzales, P., Frase, M. (1999). \emph{Highlights from TIMSS}. National Center for Education Statistics, Washington DC.
    
    TIMSS International Study Center (1996). \emph{Third International Mathematics and Science Study: Technical Report, Volume 1 Design and Development}. Boston College: Chestnut Hill MA.
}
\examples{
data(ThirdGrade)
str(ThirdGrade)
}
\keyword{datasets}
