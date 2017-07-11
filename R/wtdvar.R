
wtdvar <- function(x, w){
    n <- length(w)
    xbarw <- sum(w*x) / sum(w)
    varw <- n / (n-1) * sum(w * (x-xbarw)^2) / sum(w)
    varw
}
