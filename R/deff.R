deff <- function(w, x=NULL, y=NULL, p=NULL, strvar=NULL, clvar=NULL, Wh=NULL, nest=FALSE, type){
    if (is.null(w)) stop("w is required for all deffs")

    if (!(type %in% c("kish", "henry", "spencer", "cr")))
    stop("type must be one of kish, henry, spencer, or cr. \n")
    if (type == "kish"){
        if (is.null(w)) stop("w is required for Kish deff")
        d <- deffK(w)
    }
    if (type == "henry"){
        if (any(is.null(w)|is.null(x)|is.null(y))) stop("w, x, and y are all required for Henry deff")
        d <- deffH(w=w, y=y, x=x)
    }
    if (type == "spencer"){
        if (any(is.null(w)|is.null(y)|is.null(p))) stop("w, y, and p are all required for Spencer deff")
        d <- deffS(w=w, y=y, p=p)
    }
    if (type == "cr"){
        if (any(is.null(w)|is.null(y))) stop("w,  and y are required for Chen-Rust deff")
        d <- deffCR(w=w, y=y, strvar=strvar, clvar=clvar, Wh=Wh)
    }
    return(d)
}
