SampStop <- function (lm.obj, formula, n1.data, yvar, n2.data, p = NULL, 
                      delta = NULL, seed = NULL) 
{
  if (!inherits(lm.obj, "lm")) 
    stop("lm.obj must be lm class object from lm regression.\n")
  if (is.null(p) == TRUE) 
    stop("p (response probability) must be provided.\n")
  if (is.numeric(p) == FALSE) 
    stop("p (response probability) must be numeric.\n")
  if (any(p <= 0, p >= 1)) 
    stop("p (response probability) must be greater than 0 and less than 1.\n")
  n1 <- nrow(n1.data)
  n2 <- nrow(n2.data)
  n <- n1 + n2
  Z.wr.full <- stats::model.matrix(formula, data = n2.data)
  var.e <- stats::sigma(lm.obj)^2
  mm <- stats::model.matrix(lm.obj)
  Zr.xprod <- t(mm) %*% mm
  y1 <- n1.data[, yvar]
  e1 <- mean(y1)
  d1 <- length(p)
  d2 <- length(delta)
  d3 <- 7
  out <- array(rep(0, d1 * d2 * 7), dim = c(d1, d2, 7))
  for (k in 1:length(p)) {
    n2.sample <- sample(n2, ceiling(n2 * p[k]))
    if (is.vector(Z.wr.full)) {
      Z.wr.sam <- Z.wr.full[n2.sample]
      zbar.wr <- mean(Z.wr.sam)
    }
    else {
      Z.wr.sam <- Z.wr.full[n2.sample, ]
      if (is.vector(Z.wr.sam)) {
        zbar.wr <- mean(Z.wr.sam)
      }
      else zbar.wr <- colMeans(Z.wr.sam)
    }
    var.e1.minus.e2 <- (var.e) * (n2 * p[k]/n)^2 * (1 + t(zbar.wr) %*% 
                                                      MASS::ginv(Zr.xprod) %*% zbar.wr + 1/(n2 * p[k]))
    z <- abs(delta)/rep(sqrt(var.e1.minus.e2), length(delta))
    prob.lt.delta <- pnorm(z) - pnorm(-z)
    out[k, , ] <- cbind(rep(p[k], length(delta)), 
                        round(rep(n2 * p[k], length(delta)), 0), 
                        round(rep(e1, length(delta)), 1), 
                        round(delta, 1), 
                        round(rep(sqrt(var.e1.minus.e2), length(delta)), 2), 
                        round(z, 2), 
                        round(prob.lt.delta, 3))
  }
  result <- NULL
  for (k in 1:d1) {
    result <- rbind(result, out[k, , ])
  }
  dimnames(result)[2] <- list(c("Pr(response)", "Exp no. resps", 
                                "y1 mean", "diff in means", "se of diff", "z-score", 
                                "Pr(smaller diff)"))
  result <- as.data.frame(result)
  
  suminfo <- list(c("No. of current respondents" = n1),
                  c("No. of current nonrespondents" = n2),
                  c("Formula" = deparse(formula(lm.obj))))

  suminfo <- as.data.frame(unlist(suminfo))
  colnames(suminfo) <- c("Result")
  
  result  <- list("Input"  = suminfo,
                  "Output" = result)
  
  return(result)
}
