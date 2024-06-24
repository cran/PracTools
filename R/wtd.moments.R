wtd.moments <- function(y, w=NULL, pop.sw=TRUE){
  n <- length(y)
  if (n < 2)
    stop("Length of y must be 2 or more. \n")
    if (pop.sw) {w <- rep(1, n)}

  ybarw <- sum(w*y) / sum(w)
  m2 <- wtdvar(y, w=w)
  m3 <- sum(w*(y - ybarw)^3) / sum(w)
  m4 <- sum(w*(y - ybarw)^4) / sum(w)

  skewness <- m3 / m2^(3/2)
  kurtosis <- m4 / m2^2

  c(m2=m2, m3=m3, m4=m4, skewness=skewness, kurtosis=kurtosis)
}
