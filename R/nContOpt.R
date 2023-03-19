nContOpt <- function(X,              ## Variable used for one-draw probabilities
                     Y   = NULL,     ## Variable used for estimation, design = PPS only
                     CV0 = NULL,     ## Desired CV
                     V0  = NULL,     ## Desired Variance
                                     ## Either CV0 or V0 must be set
                     design = NULL   ## Must be either "SRS" or "PPS"
                     ) {

  ## Confirm variance or CV is set, only one is set, and it is greater than 0
  if (sum(sapply(list(CV0, V0), is.null)) == 2)
    stop("At least one of CV0 or V0 must be specified\n")
  if (sum(sapply(list(CV0, V0), is.null)) == 0)
    stop("Only one of CV0 or V0 can be specified\n")
  if(any(CV0 <= 0, V0 <= 0))
    stop("CV0 and V0 cannot be <= 0.\n")

  ## Confirm sample design is "SRS" or "PPS"
  if(design != "SRS" & design != "PPS")
    stop("Sample design must be SRS or PPS.\n")

  ## Confirm X and Y are numeric
  if(is.null(X) == TRUE)
    stop("X (variable used for 1-draw probabilities) must be provided.\n")
  if(is.numeric(X) == FALSE)
    stop("X (variable used for 1-draw probabilities) must be numeric.\n")
  if(is.null(Y) == FALSE & is.numeric(Y) == FALSE)
    stop("Y (variable used for estimation in design = PPS) must be numeric.\n")


  ## Create vector for sample size
  N      <- length(X)
  n.sam  <- rep(0, N-1)

  ## Calculate sample size for each census/sample split

  if(design == "SRS"){

    ## Calculate ybarU parameter
    XbarU  <- mean(X)
    ## Sort X variable from highest to lowest
    X.sort <- sort(X, decreasing = TRUE)
    ## Create sample size vector
    ## i is the number of certainties, beginning with 0
    for(i in 1:(N-1)){
      S2X <- var(X.sort[i:N])
      if(is.null(V0) == FALSE){
        n.sam[i] <-  (i - 1) +
          ((N-i+1)/N)^2 * S2X / (V0 + (((N-i+1)/N)^2 * S2X) / (N-i+1))
      } else if(is.null(CV0) == FALSE){
        n.sam[i] <-  (i - 1) +
                    (((N-i+1)/N)^2 * S2X / XbarU^2) /
                     (CV0^2 + (((N-i+1)/N)^2 * S2X) / ((N-i+1) * XbarU^2))
        }
      }
  } else if(design == "PPS"){
    ## Sort X and Y pairs from highest to lowest based on X
    XY.df   <- cbind.data.frame(X, Y)
    XY.sort <- XY.df[order(XY.df$X, decreasing = TRUE), ]
    ## Calculate ybarU for CV0 calculation
    ybarU  <- mean(Y)

    ## Create sample size vector
    ## i is the number of certainties, beginning with 0
    for(i in 1:(N-1)){
      pi  <- XY.sort$X[i:N]/sum(XY.sort$X[i:N])
      Yi  <- XY.sort$Y[i:N]
      tU  <- sum(Yi)
      S2Y <- (1 / (N + 1 - i)^2) * (sum(pi * (Yi/pi - tU)^2))
      if(is.null(V0) == FALSE){
        n.sam[i] <- (i - 1) +
                          ((N-i+1)/N)^2 * S2Y /
                    (V0 + ((N-i+1)/N)^2 * S2Y / (N-i+1))
      } else if(is.null(CV0) == FALSE){
        n.sam[i] <- (i - 1) +
                           (((N-i+1)/N)^2 * S2Y / ybarU^2) /
                   (CV0^2 + ((N-i+1)/N)^2 * S2Y / ((N-i+1) * ybarU^2))
      }
    }
  }

  Take.alls    <-       if(design == "SRS"){(X >  X.sort[which.min(n.sam)])
                 } else if(design == "PPS"){(X > XY.sort$X[which.min(n.sam)])}
  Min.Takeall.Val <-    if(design == "SRS"){(X.sort[which.min(n.sam) - 1])
                 } else if(design == "PPS"){(XY.sort$X[which.min(n.sam) - 1])}

  ## Output data for curve plot, identify take-alls, provide sample size,
  ## minimum sample value, and number of take-alls
  out <- list("nContOpt.Curve"  = n.sam,
              "Take.alls"       = Take.alls,
              "nContOpt.n"      = round(min(n.sam), 4),
              "Min.Takeall.Val" = Min.Takeall.Val,
              "n.Take.all"      = which.min(n.sam) - 1
              )

  ## Return output
  return(out)
}
