\name{GeoDistMOS}
\alias{GeoDistMOS}
\title{
   Split geographic PSUs based on a measure of size threshold
}
\description{
   Split geographic PSUs into new geographically contiguous PSUs based on a maximum measure of size for each PSU
}
\usage{
   GeoDistMOS(lat, long, psuID, n, MOS.var, MOS.takeall = 1, Input.ID = NULL)
}
\arguments{
   \item{lat}{latitude variable in an input file. Must be in decimal format.}
   \item{long}{longitude variable in an input file. Must be in decimal format.}
   \item{psuID}{PSU Cluster ID from an input file.}
   \item{n}{Sample size of PSUs; may be a preliminary value used in the computation to identify certainty PSUs}
   \item{MOS.var}{Variable used for probability proportional to size sampling}
   \item{MOS.takeall}{Threshold relative measure of size value for certainties; must satisfy 0 < \code{MOS.takeall} <= 1}
   \item{Input.ID}{ID variable from the input file}
}
\details{
\code{GeoDistMOS} splits geographic primary sampling units (PSUs) in the input object based on a variable which is used to create the measure of size for each PSU (\code{MOS.var}). The goal is to create PSUs of similarly sized MOS. The input file should have one row for each geographic unit, i.e. secondary sampling unit (SSU), with a PSU ID assigned. The latitude and longitude input vectors define the centroid of each input SSU. The complete linkage method for clustering is used. Accordingly, PSUs are split on a distance metric and not on the MOS threshold value. \code{GeoDistMOS} calls the function \code{inclusionprobabilities} from the \code{sampling} package to calculate the inclusion probability for each SSU within a PSU and \code{distHaversine} from the \code{geosphere} package to calculate the distances between centroids.
}
\value{
  A list with two components:
     \item{PSU.ID.Max.MOS}{A data frame containing the SSU ID value in character format (\code{Input.ID}), the original PSU ID (\code{psuID.orig}), and the new PSU ID after splitting for the maximum measure of size (\code{psuID.new}).}
     \item{PSU.Max.MOS.Info}{A data frame containing the new PSU ID (\code{psuID.new}) after splitting for the maximum Measure of Size, the inclusion probability of the PSU ID given the input sample size n (\code{psuID.prob}), the measure of size of the new PSU (\code{MOS}), the number of SSUs in the new PSU ID (\code{Number.SSUs}), and the means of the SSUs latitudes and longitudes that were combined to form the new PSU (\code{PSU.Mean.Latitude} and \code{PSU.Mean.Longitude}).}
}

\author{
George Zipf, Richard Valliant
}
\seealso{
\code{\link{GeoDistPSU}}, \code{\link{GeoMinMOS}}
}
\examples{
data(Test_Data_US)

   # Create PSU ID with GeoDistPSU
g <- GeoDistPSU(Test_Data_US$lat,
                Test_Data_US$long,
                "miles",
                100,
                Input.ID = Test_Data_US$ID)
   # Append PSU ID to input file
library(dplyr)
Test_Data_US <- dplyr::inner_join(Test_Data_US, g$PSU.ID, by=c("ID" = "Input.file.ID"))

   # Split PSUs with MOS above 0.80
m <- GeoDistMOS(lat         = Test_Data_US$lat,
                long        = Test_Data_US$long,
                psuID       = Test_Data_US$psuID,
                n           = 15,
                MOS.var     = Test_Data_US$Amount,
                MOS.takeall = 0.80,
                Input.ID    = Test_Data_US$ID)

   # Create histogram of Measure of Size Values
hist(m$PSU.Max.MOS.Info$psuID.prob,
     breaks = seq(0, 1, 0.1),
     main = "Histogram of PSU Inclusion Probabilities (Certainties = 1)",
     xlab = "Inclusion Probability",
     ylab = "Frequency")
}
\keyword{methods}
\keyword{survey}

