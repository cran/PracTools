BW2stagePPS <- function(X, pp, psuID, lonely.SSU = "mean"){
    if (any(pp >=1)){
        pick <- (1:length(pp))[pp < 1]
        warning(paste(sum(pp >= 1), " PSUs are excluded from computations because pp >= 1.\n\n"))
        keep <- (psuID %in% pick)
        X <- X[keep]
        pp <-  pp[pick]
        psuID <- psuID[keep]
    }
    
    if (!(sum(pp)==1)) {stop("pp vector must sum to 1.\n") }

    M <- length(unique(psuID))
    Ni <- table(psuID)
    cl.tots <- by(X, INDICES = psuID, FUN = sum)
    cl.vars <- by(X, INDICES = psuID, FUN = var)
    cl.vars.miss <- is.na(cl.vars)
    if (lonely.SSU == "mean"){
      cl.vars[cl.vars.miss] <- mean(cl.vars[!cl.vars.miss])
    }
    else if (lonely.SSU == "zero"){
      cl.vars[cl.vars.miss] <- 0
    }
    else {stop("Illegal value of lonely.SSU: ", lonely.SSU, "\n")}

    tU <- sum(cl.tots)
    S2U1 <- sum(pp * (cl.tots/pp - tU)^2)
    B2 <- S2U1 / tU^2

    ybarU <- mean(X)
    W2 <- sum(Ni^2 * cl.vars / pp) / tU^2
    S2U <- var(X)
    k <- (B2 + W2) / (S2U/ybarU^2)

    c("B2"=B2,
      "W2"=W2,
      "unit relvar" = S2U/ybarU^2,
      "B2+W2"=B2 + W2,
      "k"=k,
      "delta"=B2 / (B2 + W2))
}
