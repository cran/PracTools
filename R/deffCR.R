deffCR <- function(w, strvar=NULL, clvar=NULL, Wh=NULL, nest=FALSE, y){
    if (any(is.null(w) | is.null(y))) stop("w and y are required for Chen-Rust deff")
    n <- length(w)
    sig2 <- n / (n-1) * sum(w*(y - sum(w*y)/sum(w))^2) / (sum(w)-1)

    if (!is.null(strvar)){             # strata = Y
        strat <- unique(strvar)
        H <- length(strat)
        nh <- as.vector(table(strvar))
        sig2h <- deff.s <- cv2h <- vector("numeric", length=H)
        if (is.null(Wh)){           # estimate Wh if it's not fed as a parameter
            Wh <- by(w, INDICES=strvar, sum) / sum(w)
            Wh <- as.vector(Wh)
        }
        for (h in strat){
            wsub <- w[strvar == h]
            ysub <- y[strvar == h]
            cv2h[h] <- (nh[h]-1)/nh[h] * var(wsub) / mean(wsub)^2
            if (is.null(clvar)){       # strata = Y, clusters = N
                sig2h[h] <- nh[h] / (nh[h]-1) * sum(wsub*(ysub - sum(wsub*ysub)/sum(wsub))^2) / (sum(wsub)-1)
                deff.s[h] <- Wh[h]^2/ nh[h] * n * sig2h[h] / sig2
            }
        }

        if (is.null(clvar)) {
            deff.w <- 1 + cv2h
            out <- list("strata deffs" = cbind(stratum = 1:H, deff.w = deff.w, deff.c = NULL, deff.s = deff.s),
                        "overall deff" = sum(deff.w * deff.s)
                   )
        }

        if (!is.null(clvar)){           # strata = Y, clusters = Y
            sam.dsgn <- survey::svydesign(ids = ~clvar, strata = ~strvar, data = data.frame(y), weights = w, nest = nest)
            strmns <- survey::svyby(~y, by = ~strvar, design=sam.dsgn, FUN=survey::svymean, deff=TRUE)
            strdeff <- strmns$DEff.y
            nh.star <- rhoh <- vector("numeric", length=H)

            for (h in strat){
                wsub <- w[strvar == h]
                clsub <- clvar[strvar == h]
                ysub <- y[strvar == h]
                nh.star.num <- by(wsub, INDICES=clsub, sum)
                nh.star.num <- sum(nh.star.num^2)
                nh.star.den <- sum(wsub^2)
                nh.star[h] <- nh.star.num / nh.star.den
                rhoh[h] <- (strdeff[h] - (1 + cv2h[h])) / (1 + cv2h[h]*(nh.star[h]-1))
                sig2h[h] <- nh[h] / (nh[h]-1) * sum(wsub*(ysub - sum(wsub*ysub)/sum(wsub))^2) / (sum(wsub)-1)
                deff.s[h] <- Wh[h]^2/ nh[h] * n * sig2h[h] / sig2
            }
            deff.w = 1 + cv2h
            deff.c = 1 + (nh.star-1)*rhoh
            out <- list("strata components" = cbind(stratum = 1:H, deff.w = deff.w, deff.c = deff.c, deff.s = deff.s),
                        "overall deff" = sum(deff.w * deff.c * deff.s)
                   )
        }
    }

    if (is.null(strvar)){              # strata = N
        nh <- n
        if (is.null(clvar)){           # strat = N, clusters = N
            cv2h <- (nh-1)/nh * var(w) / mean(w)^2
            deff.w = 1 + cv2h
            out <- c(deff.w = deff.w)
        }
        else {                      # strat = N, clusters = Y
            nh.star.num <- by(w, INDICES=clvar, sum)
            nh.star.num <- sum(nh.star.num^2)
            nh.star.den <- sum(w^2)
            nh.star <- nh.star.num / nh.star.den
            sam.dsgn <- survey::svydesign(ids = ~clvar, strata = NULL, data = data.frame(y), weights = w)
            strmns <- survey::svymean(~y, design=sam.dsgn, deff=TRUE)
            strdeff <- survey::deff(strmns)
            cv2h <- (nh-1)/nh * var(w) / mean(w)^2
            rhoh <- (strdeff - (1 + cv2h)) / (1 + cv2h*(nh.star-1))
            deff.w = 1 + cv2h
            deff.c = 1 + (nh.star-1)*rhoh
            out <- list(cbind(deff.w = deff.w, deff.c = deff.c, deff.s = NULL),
                        "overall deff" = deff.w * deff.c
                    )
        }
    }
    return(out)
}
