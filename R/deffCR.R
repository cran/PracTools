deffCR <- function (w, strvar = NULL, clvar = NULL, clnest = NULL, Wh = NULL, y, stages) {
  if (is.null(w) == TRUE) 
    stop("w is required for Chen-Rust deff.\n")
  if (is.null(y) == TRUE) 
    stop("y is required for Chen-Rust deff.\n")
  if (any(is.na(w)) == TRUE) 
    stop("w has missing values, which are not allowed.\n")
  if (any(is.na(y)) == TRUE) 
    stop("y has missing values, which are not allowed.\n") 
    # Check scale of weights
  if (is.null(strvar)){
    if (sum(w) <= length(w)) stop("Sum of weights is less than or equal to sample size. 
                                  deffCR requires weights that are scaled for estimating population totals.")
  }
  if (!is.null(strvar) & !is.null(clvar) & is.null(clnest)){
    stop("If sample is stratified and clustered, clnest must be either TRUE or FALSE.\n")
  }

  n <- length(w)
  sig2 <- n/(n - 1) * sum(w * (y - sum(w * y)/sum(w))^2)/(sum(w) - 1)

  if (!is.null(strvar)) {
    nh <- as.vector(table(strvar))
    ## Note: unique does not sort. Use sort(unique(strvar))
    strat <- sort(unique(strvar))
    
    for (h in strat){
      if (length(stages) != length(strat)){
      stop("The parameter stages must have the same length as the number of strata. Note that stages must be in the same order as the sorted, unique values of strvar.\n ")
      }
    }
    one.stageh <- strat[stages == 1]
    
    H <- length(strat)
    chk.w <- vector(length = H)
    
    for (h in strat) {
      indx <- (1:length(strat))[strat==h]
      chk.w[indx] <- sum(w[strvar == h]) <= length(w[strvar == h])
      if (chk.w[indx]) warning("Sum of stratum weights is less than or equal to sample size in stratum ",h,". deffCR requires weights that are scaled for estimating population totals.\n") 
    }
    if (any(chk.w)) stop("Sum of stratum weights is less than or equal to sample size in at least one stratum. 
                         Rescale weights so that they are appropriate for estimating population totals.\n")

    stagesh <- rep(stages, nh)
    sig2h <- deff.s <- cv2h <- vector("numeric", length = H)
    if (is.null(Wh)) {
      Wh <- by(w, INDICES = strvar, sum)/sum(w)
      Wh <- as.vector(Wh)
    }
    for (h in strat) {
      wsub <- w[strvar == h]
      ysub <- y[strvar == h]
      ## Issue is that cv2h is vector length H
      ## Use match function: match(h, strat)
      i <- match(h, strat)
      cv2h[i] <- (nh[i] - 1)/nh[i] * var(wsub)/mean(wsub)^2
      if (is.null(clvar)) {
        sig2h[i] <- nh[i]/(nh[i] - 1) * sum(wsub * (ysub - sum(wsub * ysub)/sum(wsub))^2)/(sum(wsub) - 1)
        deff.s[i] <- Wh[i]^2/nh[i] * n * sig2h[i]/sig2
      }
    }

    if (is.null(clvar)) {
      deff.w <- 1 + cv2h
      out <- list(`strata deffs` = data.frame(stratum = strat, nh, cv2wh = cv2h,
                                              deff.w = deff.w, deff.s = deff.s), 
                  `overall deff` = sum(deff.w * deff.s))
    }

    if (!is.null(clvar)) {
      tab.hi <- table(paste(strvar,".",clvar, sep=""))
      tab.hi <- cbind(as.numeric(substr(names(tab.hi),1,1)), tab.hi)

      stop.sw <- FALSE
      for (h in strat){
        th12 <- (as.matrix(tab.hi[tab.hi[, 1] == h, ])[,2] >= 2) & (stages[h] == 1)
        th21 <- (as.matrix(tab.hi[tab.hi[, 1] == h, ])[,2] == 1) & (stages[h] == 2)
        
        if (any(th12)) {
          stop.sw <- TRUE
          warning(paste("stages = 1 in stratum", h, " but some clusters have more than 1 element. \n"))
          }
        if (any(th21)) {
          stop.sw <- TRUE
          warning(paste("stages = 2 in stratum", h, " but some clusters have only 1 element. \n"))
          }
      }
      if (stop.sw) 
        {stop("Error: Number of stages of sampling are inconsistent with counts per cluster in one or more strata.\n")}
      
      strmns <- NULL
      for (k in 1:2){
        pick <- (stagesh == k)
        strsub <- strvar[pick]
        clsub <- clvar[pick]
        wsub <- w[pick]
        ysub <- data.frame(y = y[pick])

        if (nrow(ysub) > 0){
          if (k == 1){
            sam.dsgn <- survey::svydesign(ids = ~1, strata = ~strsub, data = ysub, weights = wsub)
            if (length(unique(strsub)) == 1){
              ysub <- as.numeric(unlist(ysub))
              mn <- survey::svymean(~ysub, design=sam.dsgn, deff = TRUE)
              strmns <- c(unique(strsub), mn, survey::SE(mn), survey::deff(mn))
            }
            else {
              ysub <- as.numeric(unlist(ysub))
              strmns <- rbind(strmns, survey::svyby(~y, by = ~strsub, design = sam.dsgn, FUN = survey::svymean, deff = TRUE, survey.lonely.psu="adjust"))              
            }
         }
          
          else if (k == 2){
            sam.dsgn <- survey::svydesign(ids = ~clsub, strata = ~strsub, 
                                          data = ysub, weights = wsub, nest = clnest)
            if (length(unique(strsub)) == 1){
              ysub <- as.numeric(unlist(ysub))
              strmns <- rbind(strmns, survey::svymean(~ysub, design=sam.dsgn, deff = TRUE))
            }
            else {
              ysub <- as.numeric(unlist(ysub))
              strmns <- rbind(strmns, survey::svyby(~y, by = ~strsub, design = sam.dsgn, FUN = survey::svymean, 
                                                    deff = TRUE, survey.lonely.psu="adjust"))              
            }
          }
        }
      }     # end k in 1:2

      strmns <- strmns[order(strmns$strsub),]
      strdeff <- strmns$DEff.y
      
      nh.star <- rhoh <- vector("numeric", length = H)

      for (h in strat) {
        wsub <- w[strvar == h]
        clsub <- clvar[strvar == h, drop = TRUE]
        ysub <- y[strvar == h]
        nh.star.num <- by(wsub, INDICES = clsub, sum)
        nh.star.num <- sum(nh.star.num^2)
        nh.star.den <- sum(wsub^2)
        i <- match(h, strat)
        nh.star[i] <- nh.star.num/nh.star.den

        if (nh.star[i] == 1){warning("Avg sample size within cluster is 1 in stratum ", i, ". rhoh is Inf and deff.c is set to 1 in stratum ",i,".\n")}
        rhoh[i]    <- (strdeff[i] - (1 + cv2h[i]))/((1 + cv2h[i]) * (nh.star[i] - 1))
        sig2h[i]   <- nh[i]/(nh[i] - 1) * sum(wsub * (ysub - sum(wsub * ysub)/sum(wsub))^2)/(sum(wsub) - 1)
        deff.s[i]  <- Wh[i]^2/nh[i] * n * sig2h[i]/sig2
      }
      deff.w = 1 + cv2h
      deff.c = 1 + (nh.star - 1) * rhoh
      deff.c[one.stageh] <- 1
      out <- list(`strata components` = data.frame(stratum = strat, nh, rhoh, cv2wh = cv2h,
                                              deff.w = deff.w, deff.c = deff.c, deff.s = deff.s), 
                  `overall deff` = sum(deff.w * deff.c * deff.s))
    }   # end !is.null(clvar)
  }   # end !is.null(strvar)
  
  if (is.null(strvar)) {
    nh <- n
    if (is.null(clvar)) {
      cv2h <- (nh - 1)/nh * var(w)/mean(w)^2
      deff.w = 1 + cv2h
      out <- list(`deff components` = data.frame(n = n, cv2w = cv2h, deff.w = deff.w),
                  `overall deff = deff.w` = deff.w)
    }
    else {
      nh.star.num <- by(w, INDICES = clvar, sum)
      nh.star.num <- sum(nh.star.num^2)
      nh.star.den <- sum(w^2)
      nh.star <- nh.star.num/nh.star.den
      sam.dsgn <- survey::svydesign(ids = ~clvar, strata = NULL, 
                                    data = data.frame(y), weights = w, nest = clnest)
      strmns <- survey::svymean(~y, design = sam.dsgn, 
                                deff = TRUE)
      strdeff <- survey::deff(strmns)
      cv2h <- (nh - 1)/nh * var(w)/mean(w)^2
      rhoh <- (strdeff - (1 + cv2h))/((1 + cv2h) * (nh.star - 1))
      deff.w = 1 + cv2h
      deff.c = 1 + (nh.star - 1) * rhoh
      out <- list(`deff components` = data.frame(n = n, rho = rhoh, cv2w = cv2h, 
                        deff.w = deff.w, deff.c = deff.c), 
                  `overall deff` = deff.w * deff.c)
    }
  }     # end is.null(strvar)
  return(out)
}
