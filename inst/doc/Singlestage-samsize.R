## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(PracTools)

## -----------------------------------------------------------------------------
    ceiling(nPropMoe(moe.sw=1, e=seq(0.01,0.08,0.01), alpha=0.05, pU=0.5))

## -----------------------------------------------------------------------------
# Neyman allocation
Nh <- c(215, 65, 252, 50, 149, 144)
Sh <- c(26787207, 10645109, 6909676, 11085034, 9817762, 44553355)
strAlloc(n.tot = 100, Nh = Nh, Sh = Sh, alloc = "neyman")

# cost constrained allocation
ch <- c(1400, 200, 300, 600, 450, 1000)
strAlloc(Nh = Nh, Sh = Sh, cost = 100000, ch = ch, alloc = "totcost")

# allocation with CV target of 0.05
strAlloc(Nh = Nh, Sh = Sh, CV0 = 0.05, ch = ch, ybarU = 11664181, alloc = "totvar")

## -----------------------------------------------------------------------------
require(PracTools)
data("smho.N874")

y <- smho.N874[,"EXPTOTAL"]
x <- smho.N874[, "BEDS"]
y <- y[x>0]
x <- x[x>0]
ybarU <- mean(y)

(N <- length(x))
CV0 <- 0.15

  # calculate V1 based on pp(x) sample
pik <- x/sum(x)
T <- sum(y)
(V1 <- sum( pik*(y/pik - T)^2))

n <- V1 / (N*ybarU*CV0)^2
(n <- ceiling(n))

  # Anticipated SE for the pps sample
(cv.pps <- sqrt(V1/(N^2*n)) / ybarU)

  # sample size for an srs to produce the same SE
ceiling(nCont(CV0 = cv.pps, S2 = var(y), ybarU = ybarU, N = N))


