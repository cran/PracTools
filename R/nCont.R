
nCont <- function(CV0=NULL, V0=NULL, S2=NULL, ybarU=NULL, N=Inf, CVpop=NULL){
    n.sam <- NULL
    if (sum(sapply(list(CV0, V0, S2, ybarU, CVpop), is.null)) == 5)
        stop("No parameters specified\n")
    
    if (sum(sapply(list(CV0, V0), is.null)) != 1)
        stop("Only one of CV0 and V0 must be specified\n")
    
    if (any(N <= 0, S2 <= 0, CV0 <= 0, V0 <=0)) 
        stop("Neither N nor S2 can be <= 0\n")

    if (sum(sapply(list(S2, ybarU, CVpop), is.null)) == 0){
        cat("S2, ybarU, and CVpop all specified. CVpop ignored.\n")
    }
    
    if (sum(sapply(list(CVpop, N, CV0), is.null)) == 0)
        n.sam <- CVpop^2 / (CV0^2 + CVpop^2/N)

    if (sum(sapply(list(S2, ybarU, N, CV0), is.null)) == 0){
        CVpop <- sqrt(S2)/ybarU
        n.sam <- CVpop^2 / (CV0^2 + CVpop^2/N)  
    }   
     
    if (sum(sapply(list(S2, N, V0), is.null)) == 0)
        n.sam <- S2 / (V0 + S2/N)

    if (is.null(n.sam)) stop("Parameter combination is wrong. Check inputs.\n")
    else n.sam
}
