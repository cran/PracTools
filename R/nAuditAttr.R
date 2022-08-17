nAuditAttr <- function(TolRate = 0.05, AccDev, CL, N = 5000){
  ## Range checks
  if(N <= 0)                      stop("N must be positive.\n")
  if(TolRate <= 0 | TolRate >= 1) stop("Tolerable rate of deviation must be between 0 and 1.\n")
  if(AccDev <  0)                 stop("The acceptable number of deviations must be positive.\n")
  if(CL <= 0 | CL >= 1)           stop("Confidence level must be between 0 and 1.\n")

  ## Set CDF vector of length = population
  n.sam <- 1:N

  ## Sample Size: Hypergeometric
  K <- floor(N * TolRate)
  cdf <- phyper(q=AccDev, m=K, n=N-K, k=n.sam)
  if (any(cdf <= 1 - CL))
       {min.nhyper <- min(n.sam[cdf <= 1 - CL])}
  else {min.nhyper <- NA}

  ## Sample Size: Binomial
  cdf <- pbinom(q=AccDev, size=n.sam, prob=TolRate)
  if (any(cdf <= 1 - CL))
       {min.nbinom <- min(n.sam[cdf <= 1 - CL])}
  else {min.nbinom <- NA}

  list("Pop.Size"                   = N,
       "Tol.Dev.Rate"               = TolRate,
       "Acceptable.Errors"          = AccDev,
       "Sample.Size.Hypergeometric" = min.nhyper,
       "Sample.Size.Binomial"       = min.nbinom
    )
}
