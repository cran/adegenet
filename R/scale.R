####################
# scaleGen methods
####################
setGeneric("scaleGen", function(x,...){standardGeneric("scaleGen")})

setMethod("scaleGen", "genind", function(x, center=TRUE, scale=TRUE,
                                      method=c("sigma", "binom"), missing=c("NA","0","mean"), truenames=TRUE){

    THRES <- 1e-10
    method <- match.arg(method)
    missing <- match.arg(missing)

    ## handle "missing" arg
    if(missing %in% c("0","mean")){
        x <- na.replace(x, method=missing, quiet=TRUE)
    }
    
    ## handle specific cases
    if(scale & tolower(method)=="binom"){
        ## get allele freq
        temp <- apply(x$tab,2,mean,na.rm=TRUE)
        ## coerce sum of alleles freq to one (in case of missing data)
        temp <- tapply(temp, x$loc.fac, function(vec) return(vec/sum(vec)))
        pbar <- unlist(temp)

        scale <- sqrt(pbar*(1-pbar))
    }

    X <- x$tab
    ## handle truenames
    if(truenames){
        X <- truenames(x)
        if(is.list(X)) { X <- X$tab }
    }
    
    ## return result
    res <- scale(X, center=center, scale=scale)
    
    ## issue a warning if some variances are null
    temp <- attr(res,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        res <- res[, !temp]
        attr(res,"scaled:center") <- attr(res,"scaled:center")[!temp]
        attr(res,"scaled:scale") <- attr(res,"scaled:scale")[!temp]
    }

    return(res)
})





setMethod("scaleGen", "genpop", function(x, center=TRUE, scale=TRUE,
                                      method=c("sigma", "binom"),  missing=c("NA","0","mean"), truenames=TRUE){

    THRES <- 1e-10
    method <- match.arg(method)
    missing <- match.arg(missing)
    
    ## make allele frequencies here
    X <- makefreq(x,quiet=TRUE,missing=missing,truenames=truenames)$tab

    ## handle specific cases
    if(scale & tolower(method)=="binom"){
        ## get allele freq
        temp <- apply(X,2,mean,na.rm=TRUE)
        ## coerce sum of alleles freq to one (in case of missing data)
        temp <- tapply(temp, x$loc.fac, function(vec) return(vec/sum(vec)))
        pbar <- unlist(temp)

        scale <- sqrt(pbar*(1-pbar))
    }

    ## return result

    res <- scale(X, center=center, scale=scale)
    
    ## issue a warning if some variances are null
    temp <- attr(res,"scaled:scale") < THRES
    if(any(temp)) {
        warning("Some scaling values are null.\n Corresponding alleles are removed.")
        res <- res[, !temp]
        attr(res,"scaled:center") <- attr(res,"scaled:center")[!temp]
        attr(res,"scaled:scale") <- attr(res,"scaled:scale")[!temp]
    }

    return(res)
})
