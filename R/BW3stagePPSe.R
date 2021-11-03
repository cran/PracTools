
BW3stagePPSe <- function(dat, v, Ni, Qi, Qij, m, lonely.SSU = "mean", lonely.TSU = "mean"){
  browser()
    y <- dat[, v]
        # 3rd stage conditional wts
    wk.ij <- dat$w / dat$w2ij
        # 3rd stage sample counts
    qij <- table(dat$ssuID)
    qbbar <- mean(qij)

        # extract first row of dat for each PSU
    xx.psu <- do.call("rbind",list(by(1:nrow(dat),dat$psuID,head,1)))
    w1i.psu <- dat[xx.psu, "w1i"]
        # compute PSU 1-draw probs
        # assume that w1i = 1/(m*pp)
    pp <- 1/(m * dat[xx.psu,]$w1i)
        # extract first row of dat for each SSU
    xx.ssu <- do.call("rbind",list(by(1:nrow(dat),dat$ssuID,head,1)))
    ni <- table(dat[xx.ssu, "psuID"])
    f2i <- ni/Ni
    nbar <- mean(ni)
    w2ij.ssu <- dat[xx.ssu, "w2ij"]

        # SSU sel probs conditional on PSU selected
    w2ijC.ssu <- w2ij.ssu/dat[xx.ssu, "w1i"]
        # estimated no. of PSUs
    M.hat <- sum(w1i.psu)

    Qbbar <- sum(w2ij.ssu*Qij) / sum(w2ij.ssu)    # estimated mean no. of elements per SSU
    Qbar <- sum(w1i.psu*Qi) / M.hat               # estimated mean no. of elements per PSU

    tij <- by(wk.ij*y, dat$ssuID, sum)
    ti <- by(as.vector(w2ijC.ssu*tij), dat[xx.ssu,]$psuID, sum)
    t.pwr <- sum(dat$w * y)

    S2ai <- by(as.vector(tij), INDICES = dat[xx.ssu,]$psuID, var)
    S2ai.miss <- is.na(S2ai)
    if (lonely.SSU == "mean"){
      S2ai[S2ai.miss] <- mean(S2ai[!S2ai.miss])
    }
    else if (lonely.SSU == "zero"){
      S2ai[S2ai.miss] <- 0
    }
    else {stop("Illegal value of lonely.SSU: ", lonely.SSU, "\n")}

    S3ij <- by(y, INDICES = dat$ssuID, var)
    V3ij <- Qij * (Qij/qij - 1) * S3ij
    V3ijb <- Qij^2 * S3ij
    S2bi <- by(as.vector(V3ij), INDICES = dat[xx.ssu,]$psuID, sum)/ni
    S2bi.miss <- is.na(S2bi)
    if (lonely.SSU == "mean"){
      S2bi[S2bi.miss] <- mean(S2bi[!S2bi.miss])
    }
    else if (lonely.SSU == "zero"){
      S2bi[S2bi.miss] <- 0
    }
    else {stop("Illegal value of lonely.SSU: ", lonely.SSU, "\n")}

    sV3i <- by(as.vector(V3ij), INDICES = dat[xx.ssu,]$psuID, sum)
    sV3ib <- by(as.vector(V3ijb), INDICES = dat[xx.ssu,]$psuID, sum)
    sV3ib.miss <- is.na(sV3ib)
    if (lonely.TSU == "mean"){
      sV3ib[!sV3ib.miss] <- mean(sV3ib[!sV3ib.miss])
    }
    else if (lonely.TSU == "zero"){
      sV3ib[!sV3ib.miss] <- 0
    }
    else {stop("Illegal value of lonely.TSU: ", lonely.TSU, "\n")}

    S1a <- sum((ti/pp - t.pwr)^2)/(m-1)
    S1b <- sum(Ni^2/ni/pp^2*((1-f2i)*S2ai + f2i*S2bi))/m
#    V3i <- by(y, INDICES = dat$psuID, FUN = wtdvar, w = dat$w)
    V3i <- vector("numeric", length = length(unique(dat$psuID)))
    PSUs <- unique(dat$psuID)
    for (ind in 1:length(PSUs)) {
      pick <- dat$psuID == PSUs[ind]
      V3i[ind] <- wtdvar(y[pick], w = dat$w[pick])
    }
    V3i.miss <- is.na(V3i)
    if (lonely.SSU == "mean"){
      V3i[V3i.miss] <- mean(V3i[!V3i.miss])
    }
    else if (lonely.SSU == "zero"){
      V3i[V3i.miss] <- 0
    }
    else {stop("Illegal value of lonely.SSU: ", lonely.SSU, "\n")}

    Vtsu <- sum(Ni^2/ni^2 * w1i.psu^2 * sV3i)
    Vssu <- sum(Ni^2/ni/(m*pp)^2*(1-f2i)*(S2ai - S2bi))
    Vpsu <- (S1a - S1b)/m

    B <-  (S1a - S1b) / t.pwr^2
    W <- sum(Qi^2 * V3i / m /pp^2) / t.pwr^2

    W2 <- sum(Ni^2/m/pp^2*(S2ai - S2bi))
    W2 <- W2 / t.pwr^2
    W3 <- sum(Ni^2/ni/m/pp^2 * sV3ib)
    W3 <- W3 / t.pwr^2

    y.mn <- sum(dat$w * y) / sum(dat$w)
    V <- wtdvar(x=y, w=dat$w)
    k1 <- (B + W)/(V/y.mn^2)
    k2 <- (W2 + W3)/(V/y.mn^2)

    delta1 <- B / (B + W)
    delta2 <- W2 / (W2 + W3)

    c(Vpsu=Vpsu, Vssu=Vssu, Vtsu=Vtsu,
      B=B, W=W, k1=k1, W2=W2, W3=W3, k2=k2,
      delta1=delta1, delta2=delta2)
}
