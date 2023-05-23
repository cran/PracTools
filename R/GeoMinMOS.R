GeoMinMOS <- function(lat,              ## Latitude variable. Must be in decimal
                      long,             ## Longitude variable. Must be in decimal
                      geo.var,          ## Geographic variable ID for grouping
                      MOS.var,          ## Variable used for PPS sampling    
                      MOS.min           ## Minimum MOS value in MOS terms
                      ) { 

  ## INPUT CHECKS
  ## Confirm latitude and longitude are numeric and have no missing values
  if(is.numeric(lat)  == FALSE)
    stop("Latitude must be numeric.\n")
  if(any(is.na(lat)) == TRUE)
    stop("Latitude has missing values, which are not allowed.\n")
  if(is.numeric(long) == FALSE)
    stop("Longiitude must be numeric.\n")
  if(any(is.na(long)) == TRUE)
    stop("Longitude has missing values, which are not allowed.\n")

  ## Confirm MOS certainty threshold is numeric and positive
  if(is.numeric(MOS.min) == FALSE)
    stop("MOS certainty threshold must be numeric.\n")
  if((MOS.min  > 0) == FALSE)
    stop("MOS certainty threshold must be greater than zero.\n")

  ## FORMAT CONVERSIONS
  ## Convert geo.var to character
  geo.var  <- as.character(geo.var)
  
  ## FUNCTION
  ## Step 1: Create PSU level data set with psuID, mean.lat, mean.long, and MOS
  ## Note: MOS var is in absolute terms, it can be, e.g., HH or $ or hospital beds
  ## The minimum MOS a set variable
  Geo.MOS       <- tapply(MOS.var, geo.var, sum)
  geo.mean.lat  <- tapply(lat,     geo.var, mean)
  geo.mean.long <- tapply(long,    geo.var, mean)
  geo.df        <- cbind(Geo.MOS, geo.mean.lat, geo.mean.long)
  geo.df        <- as.data.frame(geo.df)
  geo.df$geoID  <- rownames(geo.df)

  ## Identify PSUs below minimum required MOS threshold for merging
  Below.min.MOS <- (Geo.MOS <= MOS.min)
  geo.Below.min <- names(Below.min.MOS)[Below.min.MOS]
  
  ## Create distance matrix
  geodf <- data.frame(geo.mean.long, geo.mean.lat)
  d     <- geosphere::distm(geodf, fun = geosphere::distHaversine)
  ## Distance of d is in meters. Convert to kilometers 
  ## There are 1609.344 meters per mile
  d <- d/1000

  ## CREATE VECTORS NECESSARY FOR ANALYSIS
  ## Create lists for receiving output
  Below.min.ID <- list()       ## ID of geo.var below minimum MOS
  geo.MOS.new  <- list()       ## Incremental MOS from added geo.var  
  geo.MOS.cum  <- list()       ## Cumulative MOS to meet minimum MOS  
  geo.var.new  <- list()       ## New geo.var added to MOS
  geo.var.num  <- list()       ## Counter for each new geo.var
  geo.var.dstk <- list()       ## Distance of geo.var to new geo.var
  geo.var.dstm <- list()       ## Distance of geo.var to new geo.var
  geo.lat.new  <- list()       ## Latitude for each new geo.var
  geo.lng.new  <- list()       ## Longitude for each new geo.var
  geo.MOS.out  <- list()       ## List for output

  for (i in 1:length(geo.Below.min)) {
    ## Output variable Geo.var is the ith element of geo.Below.min vector
    ## 1) Create distance vector for each geo.var below the MOS.min
    nearest <- d[which(geo.df$geoID == geo.Below.min[i]) ,]
    ## 2) Create order vector based on distance
    ord.n   <- order(nearest)
    ## 3) Note that MOS vector is geo.MOS from above
    ## Create output data frame for each geo.Below.min
      for (j in 1:length(Geo.MOS)){
        geo.MOS.new[j]  <- Geo.MOS[ord.n[j]]
        geo.var.new[j]  <- geo.df$geoID[ord.n[j]]
        geo.var.num[j]  <- j
        geo.var.dstk[j] <- nearest[ord.n[j]]
        geo.var.dstm[j] <- nearest[ord.n[j]]/1.609344
        geo.lat.new[j]  <- geo.mean.lat[ord.n[j]]
        geo.lng.new[j]  <- geo.mean.long[ord.n[j]]
      }
      geo.MOS.cum       <- cumsum(geo.MOS.new)
      m                 <- sum(geo.MOS.cum < MOS.min) + 1

      geo.MOS.out.i     <- cbind("Geo.Var"        = geo.Below.min[i], 
                                 "New.Geo.MOS"    = geo.MOS.new,
                                 "Geo.Cum.MOS"    = geo.MOS.cum, 
                                 "Geo.Var.ID"     = geo.var.new, 
                                 "Geo.Var.Num"    = geo.var.num,
                                 "Geo.Var.Kms"    = geo.var.dstk,
                                 "Geo.Var.Miles"  = geo.var.dstm,
                                 "Geo.Var.Lat"    = geo.lat.new,
                                 "Geo.Var.Long"   = geo.lng.new)

      ## Subset data frame for one geo.var above MOS.min
      geo.MOS.out.i    <- as.data.frame(geo.MOS.out.i[1:m, ])
      geo.MOS.out[[i]] <- geo.MOS.out.i
    
  }
  geo.MOS.out <- as.data.frame(do.call(rbind, geo.MOS.out))
  geo.MOS.out <- lapply(geo.MOS.out, unlist)
  geo.MOS.out <- cbind.data.frame(geo.MOS.out)

  ## Test for duplicated combinings
  warn <- "None"
  if(any(duplicated(geo.MOS.out$Geo.Var.ID) == TRUE)){
    warning("One or more PSUs is combined with a different PSU more than once. \n")
    warn <- unique(geo.MOS.out$Geo.Var.ID[duplicated(geo.MOS.out$Geo.Var.ID)])
  }

  Parm.Info  <- cbind.data.frame("Minimum.MOS"          = MOS.min,
                                 "Geo.Vars.start"       = length(unique(geo.var)),
                                 "Geo.Vars.lt.min.MOS"  = sum(Below.min.MOS))
  Input.Info <- cbind.data.frame(Geo.MOS, Below.min.MOS)
  Input.Info$Geo.Var <- rownames(Input.Info)
  Input.Info <- Input.Info[, c(3, 1, 2)]
  rownames(Input.Info) <- NULL

  ## Output, name, and return final geo MOS output
  out        <- list(Parm.Info,
                     Input.Info,
                     geo.MOS.out,
                     warn)
  names(out) <- c("Parameter.Information",
                  "Input.Information",
                  "Geo.var.MOS.output",
                  "For.Review")
  return(out)
  
}