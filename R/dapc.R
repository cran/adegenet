#######
## dapc
########
dapc <- function (x, ...) UseMethod("dapc")

#################
## dapc.data.frame
#################
dapc.data.frame <- function(x, grp, n.pca=NULL, n.da=NULL,
                            center=TRUE, scale=FALSE, var.contrib=FALSE,
                            pca.select=c("nbEig","percVar"), perc.pca=NULL, ...){

    ## FIRST CHECKS
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!require(MASS, quiet=TRUE)) stop("MASS library is required.")
    grp <- as.factor(grp)
    if(length(grp) != nrow(x)) stop("Inconsistent length for grp")
    pca.select <- match.arg(pca.select)
    if(!is.null(perc.pca) & is.null(n.pca)) pca.select <- "percVar"
    if(is.null(perc.pca) & !is.null(n.pca)) pca.select <- "nbEig"


    ## SOME GENERAL VARIABLES
    N <- nrow(x)

    ## PERFORM PCA ##
    maxRank <- min(dim(x))

    pcaX <- dudi.pca(x, center = center, scale = scale, scannf = FALSE, nf=maxRank)
    cumVar <- 100 * cumsum(pcaX$eig)/sum(pcaX$eig)

    ## select the number of retained PC for PCA
    if(is.null(n.pca) & pca.select=="nbEig"){
            plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA")
            cat("Choose the number PCs to retain (>=1): ")
            n.pca <- as.integer(readLines(n = 1))
    }

    if(is.null(perc.pca) & pca.select=="percVar"){
        plot(cumVar, xlab="Number of retained PCs", ylab="Cumulative variance (%)", main="Variance explained by PCA")
        cat("Choose the percentage of variance to retain (0-100): ")
        nperc.pca <- as.numeric(readLines(n = 1))
    }

    ## get n.pca from the % of variance to conserve
    if(!is.null(perc.pca)){
        n.pca <- min(which(cumVar >= perc.pca))
        if(n.pca<1) n.pca <- 1
    }


    ## keep relevant PCs - stored in XU
    X.rank <- sum(pcaX$eig > 1e-14)
    n.pca <- min(X.rank, n.pca)
    if(n.pca >= N) stop("number of retained PCs of PCA is greater than N")
    if(n.pca > N/3) warning("number of retained PCs of PCA may be too large (> N /3)\n results may be unstable ")

    U <- pcaX$c1[, 1:n.pca, drop=FALSE] # principal axes
    XU <- pcaX$li[, 1:n.pca, drop=FALSE] # principal components
    XU.lambda <- sum(pcaX$eig[1:n.pca])/sum(pcaX$eig) # sum of retained eigenvalues
    names(U) <- paste("PCA-pa", 1:ncol(U), sep=".")
    names(XU) <- paste("PCA-pc", 1:ncol(XU), sep=".")


     ## PERFORM DA ##
    ldaX <- lda(XU, grp)
    if(is.null(n.da)){
        barplot(ldaX$svd^2, xlab="Linear Discriminants", ylab="F-statistic", main="Discriminant analysis eigenvalues", col=heat.colors(length(levels(grp))) )
        cat("Choose the number discriminant functions to retain (>=1): ")
        n.da <- as.integer(readLines(n = 1))
    }

    n.da <- min(n.da, length(levels(grp))-1, n.pca) # can't be more than K-1 disc. func., or more than n.pca
    predX <- predict(ldaX, dimen=n.da)


    ## BUILD RESULT
    res <- list()
    res$n.pca <- n.pca
    res$n.da <- n.da
    res$tab <- XU
    res$grp <- grp
    res$var <- XU.lambda
    res$eig <- ldaX$svd^2
    res$loadings <- ldaX$scaling[, 1:n.da, drop=FALSE]
    res$ind.coord <-predX$x
    res$grp.coord <- apply(res$ind.coord, 2, tapply, grp, mean)
    res$prior <- ldaX$prior
    res$posterior <- predX$posterior
    res$assign <- predX$class
    res$call <- match.call()

    ## optional: get loadings of alleles
    if(var.contrib){
        res$var.contr <- as.matrix(U) %*% as.matrix(ldaX$scaling)
        f1 <- function(x){
            temp <- sum(x*x)
            if(temp < 1e-12) return(rep(0, length(x)))
            return(x*x / temp)
        }
        res$var.contr <- apply(res$var.contr, 2, f1)
    }

    class(res) <- "dapc"
    return(res)
} # end dapc.data.frame





#############
## dapc.matrix
#############
dapc.matrix <- function(x, ...){
    return(dapc(as.data.frame(x), ...))
}




#############
## dapc.genind
#############
dapc.genind <- function(x, pop=NULL, n.pca=NULL, n.da=NULL,
                        scale=FALSE, scale.method=c("sigma", "binom"), truenames=TRUE, all.contrib=FALSE,
                        pca.select=c("nbEig","percVar"), perc.pca=NULL, ...){

    ## FIRST CHECKS
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!require(MASS, quiet=TRUE)) stop("MASS library is required.")

    if(!is.genind(x)) stop("x must be a genind object.")

    if(is.null(pop)) {
        pop.fac <- pop(x)
    } else {
        pop.fac <- pop
    }

    if(is.null(pop.fac)) stop("x does not include pre-defined populations, and `pop' is not provided")


    ## SOME GENERAL VARIABLES
    N <- nrow(x@tab)

    ## PERFORM PCA ##
    maxRank <- min(dim(x@tab))

    X <- scaleGen(x, center = TRUE, scale = scale, method = scale.method,
                  missing = "mean", truenames = truenames)

    ## CALL DATA.FRAME METHOD ##
    res <- dapc(X, grp=pop.fac, n.pca=n.pca, n.da=n.da,
                center=FALSE, scale=FALSE, var.contrib=all.contrib,
                pca.select=pca.select, perc.pca=perc.pca)

    res$call <- match.call()

    return(res)
} # end dapc.genind






######################
# Function print.dapc
######################
print.dapc <- function(x, ...){
    cat("\t#########################################\n")
    cat("\t# Discriminant Analysis of Principal Components #\n")
    cat("\t#########################################\n")
    cat("class: ")
    cat(class(x))
    cat("\n$call: ")
    print(x$call)
    cat("\n$n.pca:", x$n.pca, "first PCs of PCA used")
    cat("\n$n.da:", x$n.da, "discriminant functions saved")
    cat("\n$var (proportion of conserved variance):", round(x$var,3))
    cat("\n\n$eig (eigenvalues): ")
    l0 <- sum(x$eig >= 0)
    cat(signif(x$eig, 4)[1:(min(5, l0))])
    if (l0 > 5)
        cat(" ...\n\n")

    ## vectors
    sumry <- array("", c(4, 3), list(1:4, c("vector", "length", "content")))
    sumry[1, ] <- c('$eig', length(x$eig),  'eigenvalues')
    sumry[2, ] <- c('$grp', length(x$grp), 'prior group assignment')
    sumry[3, ] <- c('$prior', length(x$prior), 'prior group probabilities')
    sumry[4, ] <- c('$assign', length(x$assign), 'posterior group assignment')
    class(sumry) <- "table"
    print(sumry)

    ## data.frames
    cat("\n")
    sumry <- array("", c(5, 4), list(1:5, c("data.frame", "nrow", "ncol", "content")))
    sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "retained PCs of PCA")
    sumry[2, ] <- c("$loadings", nrow(x$loadings), ncol(x$loadings), "loadings of variables")
    sumry[3, ] <- c("$ind.coord", nrow(x$ind.coord), ncol(x$ind.coord), "coordinates of individuals (principal components)")
    sumry[4, ] <- c("$grp.coord", nrow(x$grp.coord), ncol(x$grp.coord), "coordinates of groups")
    sumry[5, ] <- c("$posterior", nrow(x$posterior), ncol(x$posterior), "posterior membership probabilities")
    class(sumry) <- "table"
    print(sumry)

    cat("\nother elements: ")
    if (length(names(x)) > 13)
        cat(names(x)[14:(length(names(x)))], "\n")
    else cat("NULL\n")
}






##############
## summary.dapc
##############
summary.dapc <- function(object, ...){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")

    x <- object
    res <- list()

    ## number of dimensions
    res$n.dim <- ncol(x$loadings)
    res$n.pop <- length(levels(x$grp))

    ## assignment success
    temp <- as.character(x$grp)==as.character(x$assign)
    res$assign.prop <- mean(temp)
    res$assign.per.pop <- tapply(temp, x$grp, mean)

    ## group sizes
    res$prior.grp.size <- table(x$grp)
    res$post.grp.size <- table(x$assign)

    return(res)
} # end summary.dapc






##############
## scatter.dapc
##############
scatter.dapc <- function(x, xax=1, yax=2, col=rainbow(length(levels(x$grp))), posi="bottomleft", bg="grey", ratio=0.3, csub=1.2, ...){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    axes <- c(xax,yax)
    par(bg=bg)
    s.class(x$ind.coord[,axes], fac=x$grp, col=col, ...)
    if(ratio>0.001) {
        add.scatter.eig(x$eig, ncol(x$loadings), axes[1], axes[2], posi=posi, ratio=ratio, csub=csub)
    }
    return(invisible(match.call()))
} # end scatter.dapc






############
## assignplot
############
assignplot <- function(x, only.grp=NULL, subset=NULL, cex.lab=.75, pch=3){
    if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")
    if(!inherits(x, "dapc")) stop("x is not a dapc object")

    if(!is.null(only.grp)){
        only.grp <- as.character(only.grp)
        ori.grp <- as.character(x$grp)
        x$grp <- x$grp[only.grp==ori.grp]
        x$assign <- x$assign[only.grp==ori.grp]
        x$posterior <- x$posterior[only.grp==ori.grp, , drop=FALSE]
    } else if(!is.null(subset)){
        x$grp <- x$grp[subset]
        x$assign <- x$assign[subset]
        x$posterior <- x$posterior[subset, , drop=FALSE]
    }


    ##table.paint(x$posterior, col.lab=ori.grp, ...)
    ## symbols(x$posterior)


    ## FIND PLOT PARAMETERS
    n.grp <- ncol(x$posterior)
    n.ind <- nrow(x$posterior)
    Z <- t(x$posterior)
    Z <- Z[,ncol(Z):1,drop=FALSE ]

    image(x=1:n.grp, y=seq(.5, by=1, le=n.ind), Z, col=rev(heat.colors(100)), yaxt="n", ylab="", xaxt="n", xlab="Clusters")
    axis(side=1, at=1:n.grp,tick=FALSE, label=colnames(x$posterior))
    axis(side=2, at=seq(.5, by=1, le=n.ind), label=rev(rownames(x$posterior)), las=1, cex.axis=cex.lab)
    abline(h=1:n.ind, col="lightgrey")
    abline(v=seq(0.5, by=1, le=n.grp))
    box()

    newGrp <- colnames(x$posterior)
    x.real.coord <- rev(match(x$grp, newGrp))
    y.real.coord <- seq(.5, by=1, le=n.ind)

    points(x.real.coord, y.real.coord, col="deepskyblue2", pch=pch)

    return(invisible(match.call()))
} # end assignplot




###############
## randtest.dapc
###############
##randtest.dapc <- function(x, nperm = 999, ...){

##} # end randtest.dapc
