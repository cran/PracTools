
BW2stagePPSe <- function(Ni, ni, X, psuID, w, m, pp){
    fi <- ni / Ni
    t.pwr <- sum(w*X)
    pi.star <- m * pp

    S2i <- by(X, INDICES = psuID, FUN = var)
    Vi <- Ni * (Ni/ni - 1) * S2i
#    F2 <- sum( (1-pi.star )/pi.star^2 * Vi ) # used in v.0.7 and earlier
    F2 <- sum(Vi /pp^2) / m^2
    F1 <- by(w*X, psuID, sum)
    v1 <- m * var(F1)

    Vpsu <- v1 - F2
    Vssu <- sum(Vi / (pi.star)^2)

#    B <- (m*v1 - sum((1-pi.star)*Vi / m / pp^2)) / t.pwr^2 # used in v.0.7 and earlier
    B <-  m*Vpsu / t.pwr^2
    W <- sum(Ni^2 * S2i / m / pp^2) / t.pwr^2
    delta <- B / (B + W)

    V <- wtdvar(x=X, w=w)
    y.mn <- sum(w*X)/sum(w)
    k <- (B + W)/(V/y.mn^2)

    c(Vpsu=Vpsu, Vssu=Vssu, B=B, W=W, k=k, delta=delta)
}
