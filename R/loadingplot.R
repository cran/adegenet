##############
# loadingplot
##############
loadingplot <- function(x, threshold=quantile(x,0.75), axis=1, fac=NULL,
                        lab=names(x), cex.lab=0.7, cex.fac=1, lab.jitter=0,
                        main="Loading plot", xlab="Variables", ylab="Loadings",...){
    ## some checks
    if(is.data.frame(x) || is.matrix(x)){
        temp <- rownames(x)
        x <- x[,axis]
        names(x) <- temp
    }
    if(!is.numeric(x)) stop("x is not numeric")
    if(any(is.na(x))) stop("NA entries in x")
    if(any(x<0)) {
        warning("Some values in x are less than 0\n Using abs(x) instead, but this might not be optimal.")
        x <- abs(x)
    }
    
    ## preliminary computations
    y.min <- min(min(x),0)
    y.max <- max(max(x),0)
    y.offset <- (y.max-y.min)*0.02
    if(is.null(lab)) {lab <- 1:length(x)}
    
    if(!is.null(fac)){
        fac <- factor(fac, levels=unique(fac))
        grp.idx <- cumsum(table(fac)) + 0.5
        grp.lab.idx <- tapply(1:length(x), fac, mean)
        grp.lab <- names(grp.idx)
        grp.idx <- grp.idx[-length(grp.idx)]
    } # end fac handling
    
    ## start the plot
    plot(x, type="h", xlab=xlab, ylab=ylab,
         main=main, xaxt="n", ylim=c(y.min,y.max*1.2), ...)

    ## add groups of variables (optional)
    if(!is.null(fac)) {
        abline(v=grp.idx,lty=2) # split groups of variables
        text(x=grp.lab.idx,y=y.max*1.15, labels=grp.lab, cex=cex.fac) # annotate groups
    }

    ## annotate variables that are above the threshold
    x.ann <- which(x > threshold)
    x.ann <- jitter(x.ann,fac=lab.jitter)
    y.ann <- x[x > threshold] + y.offset
    y.ann <- jitter(y.ann,fac=lab.jitter)
    txt.ann <- lab[x > threshold]
    text(x=x.ann, y=y.ann, label=txt.ann, cex=cex.lab)

    ## indicate the threshold
    abline(h=threshold, col="grey")

    ## build the result
    res <- list(threshold=threshold,
                var.names=txt.ann,
                var.idx=which(x > threshold),
                var.values=x[x > threshold])
    return(invisible(res))
    
} # end loadingplot
