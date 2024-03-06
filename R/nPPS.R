nPPS <- function(X   = NULL, Y   = NULL, CV0 = NULL, V0  = NULL, N   = NULL, V1  = NULL, ybarU = NULL) {

  ## Confirm variance or CV is set, only one is set, and it is greater than 0
  if (sum(sapply(list(CV0, V0), is.null)) == 2)
    stop("At least one of CV0 or V0 must be specified\n")
  if (sum(sapply(list(CV0, V0), is.null)) == 0)
    stop("Only one of CV0 or V0 can be specified\n")
  if (any(CV0 <= 0, V0 <= 0))
    stop("CV0 and V0 cannot be <= 0.\n")

  if (is.null(X) & is.null(Y)){
    if(any(is.null(N), is.null(V1), is.null(ybarU)))
        stop("If X and Y are null, N, V1, and ybarU must all be provided.\n")
  }
  if (all(!is.null(X), !is.null(Y), !is.null(V1), !is.null(ybarU), !is.null(N)))
    warning("If X, Y, V1, ybarU, and N are all provided, X and Y are ignored.\n")
    
  if (any(is.null(V1), is.null(ybarU))){
     ## Confirm X and Y are numeric and have no missing values
     if(is.null(X))
       stop("X (variable used for 1-draw probabilities) must be provided.\n")
     if(!is.numeric(X))
       stop("X (variable used for 1-draw probabilities) must be numeric.\n")
     if (any(is.na(X)))
       stop("X has missing values, which are not allowed.\n")

     if(is.null(Y))
       stop("Y (variable used for estimation) must be provided.\n")
     if(!is.numeric(Y))
       stop("Y (variable used for estimation) must be numeric.\n")
     if (any(is.na(Y)))
       stop("Y has missing values, which are not allowed.\n")

     ## Assign population value
     N     <- length(X)
     ## Calculate mean y
     ybarU <- mean(Y)
     ## Calculate V1 based on pp(X) sample
     pik   <- X/sum(X)
     T <- sum(Y)
     (V1 <- sum( pik*(Y/pik - T)^2))
}

  ## Calculate sample
  if(!is.null(CV0)) {
    n <- V1 / (N*ybarU*CV0)^2
  } else if(!is.null(V0)) {
          n <- V1 / (V0*N^2)
  }

  ## Return output
  out <- list("N" = N,
              "V1" = V1,
              "ybarU" = ybarU,
              "n" = n)
  return(out)
}
