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
            Costh = NULL,
            nh.lowerbd = NULL
){
    if (ci.lev <= 0 | ci.lev >= 1) {stop("ci.lev must be in (0,1). \n")}
    if (!(side %in% c("two.sided", "one.sided.upper", "one.sided.lower"))) {stop("side must be two.sided, one.sided.upper, or one.sided.lower. \n")}
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
          if(!is.null(Costh)){warning("Cost is ignored when allocation is Neyman.\n")}
          Costh <- NULL
        }
        if (alloc %in% c("totvar", "totcost") & (is.null(Costh)) ) {
          stop("Stratum costs Costh must be specified if allocation is totvar or totcost.\n")
        }
        if ((alloc %in% c("totcost", "totvar")))
            warning("The total sample size will be determined by Edgeworth calculation of n for normality. \n The usual textbook formulas for total sample size when the total cost or variance is constrained do not apply. \n The standard formulas for the proportions allocated to strata do apply.\n")
      }
    }
    if (!is.null(nh.lowerbd)){
      if (nh.lowerbd < 1) {stop("nh.lowerbd must be an integer larger than or equal to 1.\n")}
    }
  
    if (side == "two.sided"){
        pval <- ci.lev + (1-ci.lev)/2
    } else if (side == "one.sided.upper"){
        pval <- ci.lev
    } else if (side == "one.sided.lower"){
        pval <- ci.lev
    }

    z <- stats::qnorm(p = pval, mean = 0, sd = 1, lower.tail = TRUE)
    phi.z <- stats::dnorm(z, mean = 0, sd = 1, log = FALSE)
    
#    q.calc <- function(z, hcol, Nh, g1, g2, b2){
    q.calc <- function(z, hcol, g1, g2, b2){
        # g1, g2 = g1, g2       if no strata
        #        = gbar1, gbar2 if strata
        # b2     = 1            if no strata
        #        = b2           if strata
      q1 <- (1/6) * (2*z^2 + 1) * g1
      q2 <- (- 1/18)*(z^5 + 2*z^3 - 3*z)*(g1^2) + (1/12)*(z^3 - 3*z)*g2 - (1/4)*(z^3 + z)*b2
      names(q1) <- "q1"
      names(q2) <- "q2"
      return(c(q1, q2))
    }

    if (is.null(hcol)){   # no strata
      N <- nrow(dat)
      if (pop.sw) {   # pop
        moments <- wtd.moments(y=dat[,ycol], w=NULL, pop.sw=pop.sw)
        g1 <- moments["skewness"]
        g2 <- moments["kurtosis"] - 3
      }   # end pop
      else if (!pop.sw) {  # sample
        moments <- wtd.moments(y=dat[,ycol], w=wts, pop.sw=pop.sw)
        g1 <- moments["skewness"]
        g2 <- moments["kurtosis"] - 3
      }   # end sample

      q.12 <- q.calc(z=z, hcol=hcol, g1=g1, g2=g2, b2=1)
      
      if (side != "two.sided"){   # 1-sided CI
        a <- 1
        b <- q.12[1] * phi.z / epsilon * (1*(side=="one.sided.lower") + (-1)*(side=="one.sided.upper"))
        cc <- q.12[2] * phi.z / epsilon
        n <- max(quad_roots(a,b,cc))^2
        if (all(Im(quad_roots(a,b,cc)) !=0)) {stop(paste("All solutions for n are complex in the quadratic being solved: \n",
                                        "a=",a,"b=",b,"c=",cc,"\n"))}
      }
      else if (side == "two.sided"){   # 2-sided CI
        n <- (-2)*phi.z*q.12[2] / epsilon
      }     # end 2-sided CI
      n <- ceiling(n)
    }   # end no strata
    
    else if (!is.null(hcol)){   # strata
      H <- length(unique(dat[,hcol]))
      hlist <- unique(dat[,hcol])
      if (pop.sw) {
          Nh <- table(dat[, hcol])
          Sh <- sqrt(by(data = dat[, ycol], INDICES = dat[, hcol], FUN = var))
          g1h <- vector("numeric",length = length(hlist))
          g2h <- vector("numeric",length = length(hlist))
          h1 <- 0
          for (hh in hlist) {
            h1 <- h1 + 1
            pick <- dat[,hcol] == hh
            moments <- wtd.moments(y = dat[pick,ycol], w = NULL, pop.sw=pop.sw)
            g1h[h1] <- moments["skewness"]
            g2h[h1] <- moments["kurtosis"] - 3
          }
      }       # end pop
      else if (!pop.sw){    # sample
          Nh <- as.vector(by(wts, INDICES=dat[,hcol],sum))
          Sh <- vector("numeric",length = length(hlist))
          g1h <- vector("numeric",length = length(hlist))
          g2h <- vector("numeric",length = length(hlist))
          h1 <- 0
          for (hh in hlist) {
            h1 <- h1 + 1
            pick <- dat[,hcol] == hh
            Sh[h1] <- sqrt(wtdvar(x=dat[pick,ycol], w=wts[pick]))
            moments <- wtd.moments(y = dat[pick,ycol], w = wts[pick], pop.sw=pop.sw)
            g1h[h1] <- moments["skewness"]
            g2h[h1] <- moments["kurtosis"] - 3
          }
      }       # end sample
      
      N <- sum(Nh)
      Wh <- Nh / N

      if (alloc == "prop"){ch <- Wh}
      if (alloc == "equal"){ch <- rep(1/H, H) }
      if (alloc == "neyman"){
        ch <- Wh*Sh / sum(Wh*Sh)
      }
      if ((alloc == "totcost") | (alloc == "totvar")){
        ch <- (Wh*Sh/sqrt(Costh)) / sum((Wh*Sh/sqrt(Costh)))
      }

      b1h.num <- (Nh*Sh)^3/ch^2
      b1h.denom <- (sum( (Nh*Sh)^2/ch ))^(3/2)
      b1h <- b1h.num/b1h.denom
      gbar1 <- sum(b1h * g1h)
      
      b2h.num <- (Nh*Sh)^4/ch^3
      b2h.denom <- (sum( (Nh*Sh)^2/ch ))^2
      b2h <- b2h.num/b2h.denom
      gbar2 <- sum(b2h * g2h)
      b2 <- sum(b2h)
  
      q.12 <- q.calc(z=z, hcol=hcol, g1=gbar1, g2=gbar2, b2=b2)  
 
      if (side != "two.sided"){   # 1-sided CI
        a <- 1
        b <- q.12[1] * phi.z / epsilon * (1*(side=="one.sided.lower") + (-1)*(side=="one.sided.upper"))
        cc <- q.12[2] * phi.z / epsilon
        n <- max(quad_roots(a,b,cc))^2
        if (all(Im(n) !=0)) {stop(paste("All solutions for n are complex in the quadratic being solved: \n",
                                          "a=",a,"b=",b,"c=",cc,"\n"))}
      }   # end 1-sided CI
 
      else if (side == "two.sided"){   # 2-sided CI
        n <- (-2)*phi.z*q.12[2] / epsilon
      }     # end 2-sided CI
      n <- ceiling(n)
      nh <- round(n*ch,0)
      
      if (alloc == "totcost" | alloc == "totvar"){
        strat.dat <- cbind("Nh" = Nh, "nh" = nh, "ch" = round(ch,3),g1h, g2h, Costh)
      }
      else{
        strat.dat <- cbind("Nh" = Nh, "nh" = nh, "ch" = round(ch,3),g1h, g2h)
      }

      if (any(strat.dat[, "nh"] > Nh)){
        bad <- (1:H)[strat.dat[,"nh"] > Nh]
        strat.dat[,"nh"][bad] <- Nh[bad]
        warning("Some stratum sample allocations are > Nh and were reset to Nh.\n")
      }
      
      if (any(strat.dat[, "nh"] < nh.lowerbd) | any(strat.dat[, "nh"] < 2)) {
        if (is.null(nh.lowerbd)) {
          warning("Some stratum sample allocations are < 2. You may need to combine strata before determining an allocation or set the lowerbd parameter to be a desired minimum.\n")
        }
        if (!is.null(nh.lowerbd)) {
          bad <- (1:H)[strat.dat[,"nh"] < nh.lowerbd]
          strat.dat[,"nh"][bad] <- nh.lowerbd
          warning("Some stratum sample allocations were < nh.lowerbd and were set to equal nh.lowerbd. \n")
        }
      }
      n <- sum(strat.dat[,"nh"])
      names(n) <- NULL
      if (alloc == "totcost" | alloc == "totvar"){
        TotalCost <- sum(Costh * strat.dat[,"nh"])
      }
    }   # end strata

    if (is.null(hcol)){
      names(g1) <- names(g2) <- names(n) <- NULL
      out <- list("CI type" = paste(ci.lev, side, sep=" "),
           "epsilon" = epsilon,
           "Total sample size" = n,
           "allocation" = alloc,
           "skewness" = g1,
           "excess kurtosis" = g2,
           "Edgeworth polynomials" = c(q.12[1], q.12[2]) )
    }
    else{
      if (is.null(Costh)){
        out <- list("CI type" = paste(ci.lev, side, sep=" "),
             "epsilon" = epsilon,
             "Total sample size" = n,
             "allocation" = alloc,
             "wtd skewness" = gbar1,
             "wtd excess kurtosis" = gbar2,
             "Edgeworth polynomials" = c(q.12[1], q.12[2]), 
             "Stratum values" = data.frame("stratum" = hlist, strat.dat) )
      }
      if (!is.null(Costh)){
        out <- list("CI type" = paste(ci.lev, side, sep=" "),
             "epsilon" = epsilon,
             "Total sample size" = n,
             "allocation" = alloc,
             "Total cost" = TotalCost,
             "wtd skewness" = gbar1,
             "wtd excess kurtosis" = gbar2,
             "Edgeworth polynomials" = c(q.12[1], q.12[2]),
             "Stratum values" = data.frame("stratum" = hlist, strat.dat) )
      }
    }
  return(out)  
}
