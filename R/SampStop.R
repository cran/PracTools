SampStop <- function(lm.obj,         ## lm object from regression predicting y_hat based on n1.data
                                     ## Note: coefficient names are carried through regression
                                     ##    model variables cannot be transformed in lm.obj
                     formula,        ## righthand side of formula in lm.obj; excludes dependent y; no quotes
                     n1.data,        ## data frame containing completed sample; includes y and Z
                     yvar,           ## name or number of column in n1.data containing y
                     n2.data,        ## data frame containing unused sample; includes Z only
                     p = NULL,       ## Vector of response probabilities for n2 sample; 0<p<1
                     delta = NULL,   ## vector of e1-e2 differences in estimates for samples
                     seed = NULL     ## Seed for random sample of n2 to determine z.wr
                     ) {

  if (!inherits(lm.obj,"lm"))    
    stop("lm.obj must be lm class object from lm regression.\n")
  if(is.null(p) == TRUE)
    stop("p (response probability) must be provided.\n")
  if(is.numeric(p) == FALSE)
    stop("p (response probability) must be numeric.\n")
  if(any(p <= 0, p >= 1))
    stop("p (response probability) must be greater than 0 and less than 1.\n")

  n1 <- nrow(n1.data)
  n2 <- nrow(n2.data)
  n  <- n1 + n2
  Z.wr.full <- stats::model.matrix(formula, data = n2.data)

  var.e  <- stats::sigma(lm.obj)^2
  mm     <- stats::model.matrix(lm.obj)
  Zr.xprod <- t(mm)%*%mm
  y1 <- n1.data[, yvar]
  e1 <- mean(y1)

  d1 <- length(p)
  d2 <- length(delta)
  d3 <- 7
  out <- array(rep(0, d1*d2*7), dim=c(d1, d2, 7))

  for (k in 1:length(p)){
    n2.sample  <- sample(n2, ceiling(n2*p[k]))
    if (is.vector(Z.wr.full)) {
      Z.wr.sam <- Z.wr.full[n2.sample]
      zbar.wr <- mean(Z.wr.sam)
    } else {
      Z.wr.sam <- Z.wr.full[n2.sample,]
      if (is.vector(Z.wr.sam)) {
        zbar.wr <- mean(Z.wr.sam)
      } else zbar.wr   <- colMeans(Z.wr.sam)
    }
    
    var.e1.minus.e2 <- (var.e) * (n2*p[k]/n)^2 * (1 + t(zbar.wr) %*% MASS::ginv(Zr.xprod) %*% zbar.wr  + 1/(n2*p[k]))
    z <- abs(delta) / rep(sqrt(var.e1.minus.e2), length(delta) )
    prob.lt.delta <- pnorm(z) - pnorm(-z)

    out[k,,] <- cbind(rep(p[k], length(delta)),
                      round(rep(n2*p[k], length(delta)),0),
                      round(rep(e1, length(delta)),1),
                      round(delta,1),
                      round(rep(sqrt(var.e1.minus.e2),length(delta)),2),
                      round(z,2),
                      round(prob.lt.delta, 3)
    )
}
    result <- NULL
    for (k in 1:d1){
    result <- rbind(result,out[k,,])
    }

  dimnames(result)[2] <- list(c("Pr(response)",
                  "Exp no. resps",
                  "y1 mean",
                  "diff in means",
                  "se of diff",
                  "z-score",
                  "Pr(smaller diff)")
                  )

  return(result)
}

