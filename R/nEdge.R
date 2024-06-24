nEdge <- function(
            ci.lev,
            side,
            epsilon = 0.005,
            dat,
            pop.sw = TRUE,
            wts = NULL,
            hcol = NULL,
            ycol,
            alloc = NULL,
            Ch = NULL
){
    if (ci.lev <= 0 | ci.lev >= 1) {stop("ci.lev must be in (0,1). \n")}
    if (!(side %in% c("two.sided", "one.sided"))) {stop("side must be two.sided or one.sided. \n")}
    if (epsilon <= 0 | epsilon >= 1) {stop("epsilon must be in (0,1). \n")}

    if (is.null(alloc) & !is.null(hcol)){
        stop("Alloc must be specified if hcol is non-null. \n")
    }
    if (!is.null(alloc)){
      if (is.null(hcol))
        warning("If hcol=NULL, non-null alloc is ignored.\n")
      if (!is.null(hcol)){
        if (!(alloc %in% c("prop", "equal", "neyman", "totcost", "totvar")))
          stop("Illegal allocation specified.\n")
        if (alloc=="neyman") {
          if(!is.null(Ch)){warning("Cost is ignored when allocation is Neyman.\n")}
        }
        if (alloc %in% c("totvar", "totcost") & (is.null(Ch)) ) {
          stop("Stratum costs Ch must be specified if allocation is totvar or totcost.\n")
        }
        if ((alloc %in% c("totcost", "totvar")))
            warning("The total sample size will be determined by Edgeworth calculation of n for normality. \n The usual textbook formulas for total sample size when the total cost or variance is constrained do not apply. \n The standard formulas for the proportions allocated to strata do apply.\n")
      }
    }
    if (side == "two.sided"){
        pval <- ci.lev + (1-ci.lev)/2
    } else {
        pval <- ci.lev
    }
    z <- stats::qnorm(p = pval, mean = 0, sd = 1, lower.tail = TRUE)
    phi.z <- stats::dnorm(z, mean = 0, sd = 1, log = FALSE)

    if (is.null(hcol)){
      N <- nrow(dat)
      if (pop.sw) {
          g1 <- wtd.moments(y=dat[,ycol], w=NULL, pop.sw=pop.sw)[4]
      }
      else if (!pop.sw) {
          g1 <- wtd.moments(y=dat[,ycol], w=wts, pop.sw=pop.sw)[4]
      }
      n <- ((2*z^2 + 1) * phi.z)^2 * g1^2 / (36*epsilon^2)
      n <- round(n,0)
    }
    else if (!is.null(hcol)){
      H <- length(unique(dat[,hcol]))
      hlist <- unique(dat[,hcol])
      if (pop.sw) {
          Nh <- table(dat[, hcol])
          Sh <- sqrt(by(data = dat[, ycol], INDICES = dat[, hcol], FUN = var))
          g1h <- vector("numeric",length = length(hlist))
          h1 <- 0
          for (hh in hlist) {
            h1 <- h1 + 1
            pick <- dat[,hcol] == hh
            g1h[h1] <- wtd.moments(y = dat[pick,ycol], w = NULL, pop.sw=pop.sw)[4]
          }
      }
      else if (!pop.sw){
          Nh <- as.vector(by(wts, INDICES=dat[,hcol],sum))
          Sh <- vector("numeric",length = length(hlist))
          g1h <- vector("numeric",length = length(hlist))
          h1 <- 0
          for (hh in hlist) {
            h1 <- h1 + 1
            pick <- dat[,hcol] == hh
            Sh[h1] <- sqrt(wtdvar(x=dat[pick,ycol], w=wts[pick]))
            g1h[h1] <- wtd.moments(y = dat[pick,ycol], w = wts[pick], pop.sw=pop.sw)[4]
          }
      }
      N <- sum(Nh)
      Wh <- Nh / N

      if (alloc == "prop"){ph <- Wh}
      if (alloc == "equal"){ph <- rep(1/H, H) }
      if (alloc == "neyman"){
        ph <- Wh*Sh / sum(Wh*Sh)
      }
      if ((alloc == "totcost") | (alloc == "totvar")){
        ph <- (Wh*Sh/sqrt(Ch)) / sum((Wh*Sh/sqrt(Ch)))
      }

      bh.num <- (Nh*Sh)^3/ph^2
      bh.denom <- (sum( (Nh*Sh)^2/ph ))^(3/2)
      bh <- bh.num/bh.denom

      gbar1 <- sum(bh * g1h)
      n <- ((2*z^2 + 1) * phi.z)^2 * gbar1^2 / (36*epsilon^2)
      n <- round(n,0)
      strat.dat <- cbind("nh" = round(n*ph,0), "ph" = round(ph,3),g1h)
      if (any(strat.dat[, "nh"] < 2)) warning("Some stratum sample allocations are < 2. You may need to combine strata before determining an allocation.\n")
      if (any(strat.dat[, "nh"] > Nh)) warning("Some stratum sample allocations are > Nh.\n")
    }

    if (is.null(hcol)){
      list("CI type" = paste(ci.lev, side, sep=" "),
           "epsilon" = epsilon,
           "Total sample size" = n,
           "allocation" = alloc,
           "g1" = g1)
    }
    else{
        list("CI type" = paste(ci.lev, side, sep=" "),
             "epsilon" = epsilon,
             "Total sample size" = n,
             "allocation" = alloc,
             "Stratum values" = data.frame("stratum" = hlist, strat.dat))
    }
}
