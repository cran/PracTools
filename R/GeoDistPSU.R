GeoDistPSU <- function(lat,              ## Latitude variable. Must be in decimal
                       long,             ## Longitude variable. Must be in decimal
                       dist.sw,          ## Distance: miles or kilometers (kms)
                       max.dist,         ## Maximum distance within PSU
                       Input.ID = NULL   ## ID variable from input file
                       ) {

  ## Confirm latitude and longitude are numeric and have no missing values
  if(is.numeric(lat)  == FALSE)
    stop("Latitude must be numeric.\n")
  if(any(is.na(lat)) == TRUE)
    stop("Latitude has missing values, which are not allowed.\n")

  if(is.numeric(long) == FALSE)
    stop("Longiitude must be numeric.\n")
  if(any(is.na(long)) == TRUE)
    stop("Longitude has missing values, which are not allowed.\n")

  ## Confirm distance switch is "miles" or "kms"
  if(dist.sw != "miles" & dist.sw != "kms")
    stop("Distance switch must be miles or kms (kilometers).\n")

  ## Confirm distance is numeric and positive
  if(is.numeric(max.dist) == FALSE)
    stop("Maximum distance must be numeric.\n")
  if((max.dist > 0) == FALSE)
   stop("Maximum distance must be greater than zero.\n")

  ## Create distance matrix from latitude and longitude
  geodf <- data.frame(long, lat)
  ## Create "As the crow flies" distance with lat and long
  d <- geosphere::distm(geodf, fun = geosphere::distHaversine)
  ## Distance of d is in meters
  ## There are 1609.344 meters per mile
  d <- if(dist.sw == "miles") d/1609.344 else
                              d/1000
  ## Convert d to distance matrix
  dist <- as.dist(d)
  ## Perform hierarchical clustering using maximum distance between two cluster objects
  hc   <- hclust(dist, method = "complete")
  ## Cut dendogram by maximum distance for PSU assignment
  psuID  <- cutree(hc, h = max.dist)

  ## Create plot of PSU centers
  PSU.Mean.Latitude  <- tapply(lat,  psuID, mean)
  PSU.Mean.Longitude <- tapply(long, psuID, mean)

  ## Calculate maximum distance between units within cluster
  PSU.Max.Dist <- NULL
  for(i in 1:length(unique(psuID))){
    PSU.Max.Dist[i] <- max(geosphere::distm(geodf[psuID == i, ],
                                                fun = geosphere::distHaversine))}
  PSU.Max.Dist <- if(dist.sw == "miles") PSU.Max.Dist/1609.344 else
    PSU.Max.Dist/1000

  ## Calculate number of SSUs in each PSU
  Number.SSUs <- table(psuID)

  ## Carry Input.ID through
  Input.file.ID <- if(is.null(Input.ID)) seq(1, length(psuID)) else
    Input.ID

  ## Create data frame with Input file ID and PSU Cluster ID
  PSU.ID   <- cbind(Input.file.ID,
                    psuID)
  PSU.ID   <- as.data.frame(PSU.ID)

  ## Create data frame with PSU Centroids, Number of SSUs, and Maximum Cluster Distance
  PSU.Info <- cbind(Number.SSUs,
                    PSU.Mean.Latitude,
                    PSU.Mean.Longitude,
                    PSU.Max.Dist)
  PSU.Info <- as.data.frame(PSU.Info)

  ## Output PSU ID and PSU Information data frames
  out <- list(PSU.ID,
              PSU.Info)

  ## Name data frames in list
  names(out) <- c("PSU.ID", "PSU.Info")

  ## Return output
  return(out)
  }
