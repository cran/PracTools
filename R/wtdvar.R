wtdvar <- function(x, w, na.rm = TRUE){
    if (na.rm){
        isNA <- is.na(x) | is.na(w)
        x <- x[!isNA]
        w <- w[!isNA]
    }
    n <- length(w)
    xbarw <- sum(w*x) / sum(w)
    varw <- n / (n-1) * sum(w * (x-xbarw)^2) / sum(w)
    varw
}
