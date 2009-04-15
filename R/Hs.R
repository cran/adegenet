############################
# Hs (expected heterozygosity)
############################
Hs <- function(x, truenames=TRUE) {

    ## checks
    if(!is.genpop(x)) stop("x is not a valid genpop object")
    if(x@type=="PA") stop("not implemented for presence/absence markers")

    ## main computations
    x.byloc <- seploc(x, truenames=truenames)
    lX <- lapply(x.byloc, function(e) makefreq(e, quiet=TRUE, truenames=truenames)$tab)
    lres <- lapply(lX, function(X) 1- apply(X^2,1,sum))

    res <- apply(as.matrix(data.frame(lres)),1,mean)

    return(res)
} # end Hs
