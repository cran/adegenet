#################
# fstat function
#################
#
# Wrapper for fst estimator from hierfstat package
#
fstat <- function(x, pop=NULL, fstonly=FALSE){
    ## misc checks
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")
    if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
    checkType(x)

    if(is.null(pop)) pop <- x@pop
    if(is.null(pop)) stop("no pop factor provided")
    if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")

    ## computations
    dat <- genind2hierfstat(x)[,-1]
    res <- varcomp.glob(levels=data.frame(pop), loci=dat)$F

    if(fstonly) {res <- res[1,1]}
    return(res)
}



## ###############
## # fst function
## ###############
## #
## # classical fst sensu Weir 1996 Genetic data analysis II pp. 166-167
## #
## fst <- function(x, pop=NULL){
##     ## misc checks
##     if(!is.genind(x)) stop("x is not a valid genind object")
##     if(!require(hierfstat)) stop("hierfstat package is required. Please install it.")

##     if(is.null(pop)) pop <- x@pop
##     if(is.null(pop)) stop("no pop factor provided")
##     if(length(pop)!=nrow(x@tab)) stop("pop has a wrong length.")

##     ## computations

##     return(res)
## }
