#################
# fstat function
#################
fstat <- function(x, pop=NULL, fstonly=FALSE){
    ## misc checks
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")
    
    if(is.null(pop)) pop <- x@pop
    if(is.null(pop)) stop("no pop factor provided")
    if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")

    ## computations
    dat <- genind2hierfstat(x)[,-1]
    res <- varcomp.glob(levels=data.frame(pop), loci=dat)$F

    if(fstonly) {res <- res[1,1]}
    return(res)
}
