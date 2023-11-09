
deffCR <- function (w, strvar = NULL, clvar = NULL, Wh = NULL, nest = FALSE, y) {

  if (is.null(w) == TRUE) 
    stop("w is required for Chen-Rust deff.\n")
  if (is.null(y) == TRUE) 
    stop("y is required for Chen-Rust deff.\n")
  if (any(is.na(w)) == TRUE) 
    stop("w has missing values, which are not allowed.\n")
  if (any(is.na(y)) == TRUE) 
    stop("y has missing values, which are not allowed.\n") 

  n <- length(w)
  sig2 <- n/(n - 1) * sum(w * (y - sum(w * y)/sum(w))^2)/(sum(w) - 1)
  if (!is.null(strvar)) {
    ## Note: unique does not sort. Use sort(unique(strvar))
    strat <- sort(unique(strvar))
    H <- length(strat)
    ## Note: table does sort
    nh <- as.vector(table(strvar))
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
                                   deff.w = deff.w, deff.c = NULL, deff.s = deff.s), 
                  `overall deff` = sum(deff.w * deff.s))
    }
    if (!is.null(clvar)) {
      sam.dsgn <- survey::svydesign(ids = ~clvar, strata = ~strvar, 
                                    data = data.frame(y), weights = w, nest = nest)
      strmns <- survey::svyby(~y, by = ~strvar, design = sam.dsgn, 
                              FUN = survey::svymean, deff = TRUE)
      strdeff <- strmns$DEff.y
      nh.star <- rhoh <- vector("numeric", length = H)
      for (h in strat) {
        wsub <- w[strvar == h]
        clsub <- clvar[strvar == h]
        ysub <- y[strvar == h]
        nh.star.num <- by(wsub, INDICES = clsub, sum)
        nh.star.num <- sum(nh.star.num^2)
        nh.star.den <- sum(wsub^2)
        ## Issue is that cv2h is vector length H
        ## Use match function: match(h, strat)
        i <- match(h, strat)
        nh.star[i] <- nh.star.num/nh.star.den
        rhoh[i]    <- (strdeff[i] - (1 + cv2h[i]))/((1 + cv2h[i]) * (nh.star[i] - 1))
        sig2h[i]   <- nh[i]/(nh[i] - 1) * sum(wsub * (ysub - sum(wsub * ysub)/sum(wsub))^2)/(sum(wsub) - 1)
        deff.s[i]  <- Wh[i]^2/nh[i] * n * sig2h[i]/sig2
      }
      deff.w = 1 + cv2h
      deff.c = 1 + (nh.star - 1) * rhoh
      out <- list(`strata components` = data.frame(stratum = strat, nh, rhoh, cv2wh = cv2h,
                                              deff.w = deff.w, deff.c = deff.c, deff.s = deff.s), 
                  `overall deff` = sum(deff.w * deff.c * deff.s))
    }
  }
  if (is.null(strvar)) {
    nh <- n
    if (is.null(clvar)) {
      cv2h <- (nh - 1)/nh * var(w)/mean(w)^2
      deff.w = 1 + cv2h
      out <- data.frame(n = n, cv2w = cv2h, deff.w = deff.w)
    }
    else {
      nh.star.num <- by(w, INDICES = clvar, sum)
      nh.star.num <- sum(nh.star.num^2)
      nh.star.den <- sum(w^2)
      nh.star <- nh.star.num/nh.star.den
      sam.dsgn <- survey::svydesign(ids = ~clvar, strata = NULL, 
                                    data = data.frame(y), weights = w)
      strmns <- survey::svymean(~y, design = sam.dsgn, 
                                deff = TRUE)
      strdeff <- survey::deff(strmns)
      cv2h <- (nh - 1)/nh * var(w)/mean(w)^2
      rhoh <- (strdeff - (1 + cv2h))/((1 + cv2h) * (nh.star - 1))
      deff.w = 1 + cv2h
      deff.c = 1 + (nh.star - 1) * rhoh
      out <- list(data.frame(n = n, rho = rhoh, cv2w = cv2h, 
                        deff.w = deff.w, deff.c = deff.c, 
                        deff.s = NULL), 
                  `overall deff` = deff.w * deff.c)
    }
  }
  return(out)
}

