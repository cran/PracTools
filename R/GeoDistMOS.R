
GeoDistMOS <- function(lat,              ## Latitude variable. Must be in decimal
                       long,             ## Longitude variable. Must be in decimal
                       psuID,            ## PSU Cluster ID
                       n,                ## Sample size
                       MOS.var,          ## Variable used for PPS sampling    
                       MOS.takeall = 1,  ## Threshold value for MOS certainties
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

  ## Confirm sample size is numeric and positive 
  if(is.numeric(n) == FALSE)
    stop("Sample size must be numeric.\n")
  if((n > 0) == FALSE)
   stop("Sample size must be greater than zero.\n")

  ## Confirm MOS certainty threshold is numeric, positive, and less than 1
  if(is.numeric(MOS.takeall) == FALSE)
    stop("MOS certainty threshold must be numeric.\n")
  if((MOS.takeall  > 0) == FALSE)
    stop("MOS certainty threshold must be greater than zero.\n")
  if((MOS.takeall <= 1) == FALSE)
    stop("MOS certainty threshold must be less than or equal to one.\n")
  
  ## Convert psuID to character
  psuID    <- as.character(psuID)
  ## Convert Input.ID to character
  Input.ID <- as.character(Input.ID)
  
  ## Step 1: Identify certainty PSUs
  ## Calculate MOS for each PSU from MOS variable
  psuID.MOS  <- (tapply(MOS.var, psuID, sum) * n) / sum(MOS.var)
  ## Compute inclusion probabilities based on PSU MOS and add PSU name
  psuID.pik  <- sampling::inclusionprobabilities(psuID.MOS, n)
  names(psuID.pik) <- names(psuID.MOS)
  ## Identify PSUs above MOS certainty threshold for splitting
  Cert       <- (psuID.pik >= MOS.takeall)
  psuID.cert <- names(Cert)[Cert]
  
  ## Step 2: Identify certainty sampling units within certainty PSU
  ## Note: psuID.pik/psuID.MOS gives adjustment factor from MOS to pik
  psu.MOS     <- list()
  psu.MOS.adj <- psuID.pik/psuID.MOS
  psu.pik     <- list()
  psu.cert    <- list()
  for (i in psuID.cert){
    psu.MOS[[i]]   <- MOS.var[psuID == i] * n / sum(MOS.var)
    psu.pik[[i]]   <- psu.MOS[[i]] * psu.MOS.adj[i]
    psu.cert[[i]]  <- (psu.pik[[i]] > MOS.takeall)
  }
  
  ## Step 3: Create hierarchical cluster based on distance for within PSU
  ## See notes for GeoDistPSU
  psu.split <- list()
  cert.ID   <- list()
  for (i in psuID.cert){
    
    ## Select psuID sampling units and create within psuID clusters
    geodf        <- data.frame(long, lat)[psuID == i, ]
    cert.ID[[i]] <- Input.ID[as.numeric(rownames(geodf))]
    ## Note: geodf has the row numbers of the original data frame. 
    ## Input ID can be compared with that
    d       <- geosphere::distm(geodf, fun = geosphere::distHaversine)
    dist    <- as.dist(d)
    hc      <- hclust(dist, method = "complete")

    ## Create psu.cert for each certainty sampling unit
    psu.cert[[i]]  <- psu.cert[[i]] * cumsum(psu.cert[[i]])
    ## Create psuID.split for clustering
    psu.split[[i]] <- psu.cert[[i]]
    ## If the sum of the non-certainties are less than MOS.takeall, put in separate psuID
    if(sum(cutree(hc, k=1)[psu.cert[[i]] == 0] * 
           psu.pik[[i]][psu.cert[[i]] == 0]) < MOS.takeall) {
      psu.split[[i]][psu.cert[[i]] == 0] <- max(psu.cert[[i]]) + 1
    } 
    else {
      ## Identify optimum k for MOS.takeall where non-certainties need to be split
      k <- 1
      while(max(tapply(psu.pik[[i]][psu.cert[[i]] == 0], 
                       cutree(hc, k = k)[psu.cert[[i]] == 0], 
                       sum)) >= MOS.takeall){
        k <- k + 1 
        if(max(tapply(psu.pik[[i]][psu.split[[i]] == 0], 
                      cutree(hc, k = k)[psu.split[[i]] == 0], 
                      sum)) < MOS.takeall){
          psu.split[[i]][psu.split[[i]] == 0] <- cutree(hc, k = k) + max(psu.split[[i]])
          }
        }
      
       }
    }

  ## Create variables for psuID, psuID.split
  ## Use sprintf type thing over psuID.cert, or length(psu.split)
  ps <- list()
  for (i in psuID.cert){
    ps[[i]]$psuID         <- psuID[psuID == i]
    ps[[i]]$psuID.max.MOS <- psu.split[which(i == psuID.cert)]
    ps[[i]]$Input.ID      <- cert.ID[i]

  }
  
  for (i in 1:length(psuID.cert)){
    ps[[i]] <- as.data.frame(ps[[i]])
    colnames(ps[[i]]) <- c("psuID.orig", "psuID.new", "Input.ID")
  }
  
  ## Create data frame of split PSU observations
  ps  <- do.call(rbind.data.frame, ps)
  ps$psuID.new  <- paste0(ps$psuID.orig, ".", ps$psuID.new)
  row.names(ps) <- NULL
  
  ## Merge with un-split PSUs
  orig.psu <- cbind.data.frame(Input.ID, psuID, MOS.var)
  ps       <- dplyr::full_join(ps, orig.psu, by=c("Input.ID"))
  ## Add psu.ID.orig and psu.ID.new to ps data frame
  ps$psuID.orig <- ifelse(is.na(ps$psuID.orig), ps$psuID, ps$psuID.orig)
  ps$psuID.new  <- ifelse(is.na(ps$psuID.new), paste0(ps$psuID.orig, ".", "0"), ps$psuID.new)
  
  ## Create MOS data frame for histogram
  mos            <- (tapply(ps$MOS.var, ps$psuID.new, sum) * n) / sum(MOS.var)
  mos            <- as.data.frame(mos) 
  mos$psuID.new  <- rownames(mos)
  rownames(mos)  <- NULL
  mos$psuID.prob <- sampling::inclusionprobabilities(mos$mos, n)
  ## Calculate number of SSUs in each PSU
  mos$Number.SSUs <- table(ps$psuID.new)
  
  ## Output data frame of PSUs split, split psuIDs, and matching information
  ## Remove un-needed variables and re-order columns
  ps  <- ps[, c(3, 1, 2)]
  mos <- mos[, c(2, 3, 4)]
  
  ## Output final MOS distribution
  out <- list(ps,
              mos)
  
  ## Name data frames in list
  names(out) <- c("PSU.ID.Max.MOS", "PSU.Max.MOS.Info")
  
  ## Return output and class
  return(out)
  
  }



