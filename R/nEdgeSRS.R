nEdgeSRS <- function(
            ci.lev,             
            side,               
            epsilon = 0.005,    
            dat,                
            pop.sw = TRUE,
            wts = NULL,
            hcol = NULL,               
            ycol                
){
    if (ci.lev <= 0 | ci.lev >= 1) {stop("ci.lev must be in (0,1). \n")}
    if (!(side %in% c("two.sided", "one.sided"))) {stop("side must be two.sided or one.sided. \n")}
    if (epsilon <= 0 | epsilon >= 1) {stop("epsilon must be in (0,1). \n")}
  
    if (side == "two.sided"){
        pval <- ci.lev + (1-ci.lev)/2
    }
    else {
        pval <- ci.lev
    }
    u <- stats::qnorm(p = pval, mean = 0, sd = 1, lower.tail = TRUE)
    phi.u <- stats::dnorm(u, mean = 0, sd = 1, log = FALSE)
    H3 <- u^3 - 3*u
    H5 <- u^5 - 10*u^3 + 15*u

    if (is.null(hcol)){
      if (pop.sw){
            N <- nrow(dat)
            g1 <- wtd.moments(y = dat[,ycol], w=NULL, pop.sw=pop.sw)[4:5]     
            g1[2] <- g1[2] - 3   # adjust to use Sugden, Smith, Jones defn
      }
      if (!pop.sw){
          N <- sum(wts)  
          g1 <- wtd.moments(y = dat[,ycol], w = wts, pop.sw=pop.sw)[4:5]
          g1[2] <- g1[2] - 3   # adjust to use Sugden, Smith, Jones defn
      }
      a <- 36*epsilon / (N*phi.u) - g1[2]*(9*H3 + 36*u)/N^2 + (g1[1]^2)*(72*u + 24*H3 + H5)/N^2
      b <- (-1)*(36*epsilon / phi.u - g1[2]*(18*H3 + 36*u)/N + (g1[1]^2)*(144*u + 72*H3 +4*H5)/N + (72*u + 18*H3)/N )
      cc <- 72*u + 18*H3 - 6*g1[2]*H3 + (g1[1]^2)*(72*u + 48*H3 + 4*H5)
      
      roots <- quad_roots(a,b,cc)
      n <- ceiling(min(roots))
      n.cochran <- ceiling(28 + 25*g1[1]^2)
      names(n.cochran) <- NULL
    }
    else if (!is.null(hcol)){
        hlist <- unique(dat[,hcol])
        if (pop.sw){ 
            Nh <- table(dat[, hcol])
            g1h <- matrix(nrow = length(hlist), ncol = 2)
            h1 <- 0
            for (hh in hlist) {
                h1 <- h1 + 1
                pick <- dat[,hcol] == hh
                g1h[h1,] <- wtd.moments(y = dat[pick,ycol], w = NULL, pop.sw=pop.sw)[4:5]
                g1h[h1,2] <- g1h[h1,2] - 3   # adjust to use Sugden, Smith, Jones defn                
            }
        }  
        else if (!pop.sw){
            Nh <- as.vector(by(wts, INDICES=dat[,hcol],sum))
            g1h <- matrix(nrow = length(hlist), ncol = 2)
            h1 <- 0
            for (hh in hlist) {
                h1 <- h1 + 1
                pick <- dat[,hcol] == hh
                g1h[h1,] <- wtd.moments(y = dat[pick,ycol], w = wts[pick], pop.sw=pop.sw)[4:5]
                g1h[h1,2] <- g1h[h1,2] - 3   # adjust to use Sugden, Smith, Jones defn
            }
        }    
        a <- 36*epsilon / (Nh*phi.u) - g1h[,2]*(9*H3 + 36*u)/Nh^2 + (g1h[,1]^2)*(72*u + 24*H3 + H5)/Nh^2
        b <- (-1)*(36*epsilon / phi.u - g1h[,2]*(18*H3 + 36*u)/Nh + (g1h[,1]^2)*(144*u + 72*H3 +4*H5)/Nh + (72*u + 18*H3)/Nh )
        cc <- 72*u + 18*H3 - 6*g1h[,2]*H3 + (g1h[,1]^2)*(72*u + 48*H3 + 4*H5)
    
        roots <- quad_roots(a,b,cc)
        roots <- matrix(roots, ncol=2, byrow = FALSE)
        roots[roots<0] <- Inf

        nh <- ceiling(apply(roots, 1, FUN = min))
        
        nh.cochran <- ceiling(28 + 25*g1h[,1]^2)

        strat.dat <- cbind(nh, round(nh/sum(nh),3), nh.cochran, g1h)
        colnames(strat.dat) <- c("nh", "proportion.nh", "nh.cochran", "stratum.skewness","stratum.kurtosis")
    }
    
    if (is.null(hcol)){
      list("CI type" = paste(ci.lev, side, sep=" "),
           "epsilon" = epsilon,
           "Total sample size" = c("n.quad"=n, "n.cochran"=n.cochran),
           "g1" = g1)
    }
    else {    
      list("Total sample size" = c("n.quad" = sum(nh), "n.cochran" = sum(nh.cochran)),
           "CI type" = paste(ci.lev, side, sep=" "),
           "epsilon" = epsilon,
           "Stratum values" = data.frame("stratum" = hlist, strat.dat))
    }
}
