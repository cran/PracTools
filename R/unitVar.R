unitVar <- function(pop.sw = NULL, w = NULL, p = NULL, y = NULL){
    if (is.null(y) & !is.numeric(y))
        stop("y must be non-null and numeric.\n")       
    if (!pop.sw & is.null(w))
        stop("If pop.sw = FALSE, w must be provided.\n")
    if (pop.sw){ 
        if (!is.null(w)) 
            warning("If pop.sw = TRUE, w is ignored.\n")
        if (!is.null(p) & length(y) != length(p))
            stop("If pop.sw = TRUE, the lengths of y and p must be equal.")
    }
    if (!is.null(w) & !is.numeric(w))
        stop("If w is provided, it must be numeric.\n")
  
    if (pop.sw & !is.null(p) & (length(y) != length(p)))
        stop ("If pop.sw = TRUE and p is provided, y is assumed to be from a sample. The lengths of p and y must be equal.\n")
    
    if (pop.sw){
        N <- length(y)
        S2 <- var(y) 
        if (!is.null(p))
            V1 <- sum( p*(y/p - sum(y))^2 ) 
    }

    if (!pop.sw) {
        S2 <- wtdvar(y, w)
        n <- length(y)
        pik <- 1/(n*w)
        if (!is.null(p)) V1 <- sum( (y/pik - sum(y/pik)/n)^2 ) / (n-1)
    }

    if (pop.sw & !is.null(p)) {
        out <- list("Note: parameters computed from full population data",
                    "Pop size N" = N,
                    "S2" = S2,
                    "V1" = V1)
    }
    if (pop.sw & is.null(p)) {out <- list("Note: parameters computed from full population data",
                      "Pop size N" = N,
                      "S2" = S2)
    }
    if (!pop.sw & !is.null(p)) out <- list("Note: parameters estimated from sample data",
                     "Sample size n" = n,
                     "S2" = S2,
                     "V1" = V1)
    if (!pop.sw & is.null(p)) out <- list("Note: parameters estimated from sample data",
                    "Sample size n" = n,
                    "S2" = S2)

    return(out)
}
