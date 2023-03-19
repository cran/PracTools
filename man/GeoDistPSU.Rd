\name{GeoDistPSU}
\alias{GeoDistPSU}
\title{
Form PSUs based on geographic distances
}
\description{
Combine geographic areas into primary sampling units to limit travel distances
}
\usage{
GeoDistPSU(lat, long, dist.sw, max.dist, Input.ID = NULL)
}
\arguments{
  \item{lat}{latitude variable in an input file. Must be in decimal format.}
  \item{long}{longitude variable in an input file. Must be in decimal format.}
  \item{dist.sw}{units for distance; either \code{"miles"} or \code{"kms"} (for kilometers)}
  \item{max.dist}{maximum distance allowed within a PSU between centroids of geographic units}
  \item{Input.ID}{ID field in the input file if present}
}
\details{
\code{GeoDistPSU} combines geographic secondary sampling units (SSUs), like cities or census block groups, into primary sampling units (PSUs) given a maximum distance allowed between the centroids of the SSUs within each grouped PSU. The input file must have one row for each geographic unit. If the input file does not have an ID field, the function will create a sequential ID that is appended to the output. The latitude and longitude input vectors define the centroid of each input SSU. The complete linkage method for clustering is used. \code{GeoDistPSU} calls the functions \code{distm} and \code{distHaversine} from the \code{geosphere} package to calculate the distances between centroids.
}
\value{
A list with two components:
    \item{PSU.ID}{A data frame with the same number of rows as the input file. Column names are \code{Input.file.ID} and \code{psuID}. The \code{psuID column} contains the PSU number assigned to each geographic unit in the input file; multiple rows of the input file will typically be assigned to the same PSU.}
    \item{PSU.Info}{A data frame with the number of rows equal to the number of PSUs that are created. Column names are \code{Num.SSUs}, number of SSUs assigned to each PSU; \code{PSU.Mean.Latitude}, mean of the latitudes of the units assigned to a PSU; \code{PSU.Mean.Longitude}, mean of the longitudes of the units assigned to a PSU; \code{PSU.Max.Dist}, maximum distance among the SSUs in a PSU}.
}

\author{
George Zipf, Richard Valliant
}
\seealso{
\code{\link{GeoDistMOS}}
}
\examples{
data(Test_Data_US)
g <- GeoDistPSU(Test_Data_US$lat,
                Test_Data_US$long,
                "miles", 100,
                Input.ID = Test_Data_US$ID)
    # Plot GeoDistPSU output
plot(g$PSU.Info$PSU.Mean.Longitude,
     g$PSU.Info$PSU.Mean.Latitude,
     col  = 1:nrow(g$PSU.Info),
     pch  = 19,
     main = "Plot of PSU Centers",
     xlab = "Longitude",
     ylab = "Latitude")
grid(col = "grey40")

    # Plot GeoDistPSU output with map
\dontrun{    
library(ggplot2)
library(sp)
library(usmap)
    # Transform PSUs into usmap projection
g.map  <- cbind(long = g$PSU.Info$PSU.Mean.Longitude,
                lat  = g$PSU.Info$PSU.Mean.Latitude)
g.map  <- as.data.frame(g.map)
g.proj <- usmap::usmap_transform(g.map,
                          input_names  = c("long", "lat"),
                          output_names = c("Long", "Lat"))
usmap::plot_usmap(color = "gray") +
  geom_point(data = g.proj,
             aes(x = Long,
                 y = Lat))
    # Create histogram of maximum distance
hist(g$PSU.Info$PSU.Max.Dist,
     main = "Histogram of Maximum Within-PSU Distance",
     xlab = "Distance",
     ylab = "Frequency")
}
}
\keyword{methods}
\keyword{survey} 