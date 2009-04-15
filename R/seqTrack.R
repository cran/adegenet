############
## generics
############
seqTrack <- function(...){
    UseMethod("seqTrack")
}

## seqTrackG <- function(...){
##     UseMethod("seqTrackG")
## }


## optimize.seqTrack <- function(...){
##     UseMethod("optimize.seqTrack")
## }


get.likelihood <- function(...){
    UseMethod("get.likelihood")
}

get.likelihood.seqTrack.default <- function(...){
    cat("Method not implemented.")
    return()
}






#######################
## auxiliary functions
#######################
#############
## .dTimeSeq
#############
##
## mu0 and L are vectors, having one value per segment/chromosome
## mu0 is per nucleotide and per day
.dTimeSeq <- function(mu, L, maxNbDays=100){
    ##mu <- mu/365 # mutation rates / site / day
    t <- 0:maxNbDays # in days added / substracted
    temp <- sapply((1-mu)^L, function(x) x^t  )
    Pt <- apply(temp,1,prod)
    t <- c(-rev(t[-1]), t)
    Pt <- c(rev(Pt[-1]), Pt)
    return(list(t, Pt))
}


#############
## .rTimeSeq
#############
##
## mu and L are vectors, having one value per segment/chromosome
##
## this returns nb days
.rTimeSeq <- function(n, mu, L, maxNbDays=100){
    temp <- .dTimeSeq(mu, L, maxNbDays)
    res <- sample(temp[[1]], size=n, replace=TRUE, prob= temp[[2]]/sum(temp[[2]]))
    return(res)
}



#################
## .rUnifDate
#################
##
## this returns random uniform dates in a given range
##
.rUnifDate <- function(n, dateMin, dateMax, ...){
    rangeSize <-  as.integer(difftime(dateMax,dateMin, units="days"))
    nbDays <- round(runif(n, min=0, max=rangeSize))
    res <- dateMin + nbDays*3600*24
    res <- as.POSIXct(round.POSIXt(res, units="days"))
    return(res)
}



#################
## .rNormTimeSeq
#################
##
## this returns nb of days
.rNormTimeSeq <- function(n, mean, sd, ...){
    res <- round(rnorm(n, mean=mean, sd=sd))
    return(res)
}



#################
## .rSampTimeSeq
#################
##
## this returns nb of days
.rSampTime <- function(n,...){
    res <- round(rnorm(n*2, -2))
    res <- res[res < 0 & res > -7][1:n]
    return(res)
}






##############
## .pAbeforeB
##############
##
## allows for different distributions for both haplo
.pAbeforeB <- function(dateA, dateB, muA, muB, LA, LB, maxNbDays=100){
    ## proba dist for A
    tempA <- .dTimeSeq(muA, LA, maxNbDays)
    days <- tempA[[1]]
    pA <- tempA[[2]]/sum(tempA[[2]]) # scale to proba mass function

    ## proba dist for B
    tempB <- .dTimeSeq(muB, LB, maxNbDays)
    pB <- tempB[[2]]/sum(tempB[[2]]) # scale to proba mass function

    ## days for A and B
    nbDaysDiff <- as.integer(round(difftime(dateA,dateB,units="days"))) # dateA - dateB, in days
    daysA <- days
    daysB <- days - nbDaysDiff

    f1 <- function(i){ # proba A before B for one day
        idx <- daysB > daysA[i]
        return(pA[i] * sum(pB[idx]))
    }

    res <- sapply(1:length(days), f1) # proba for all days
    res <- sum(res) # sum
    return(res)
}

.pAbeforeB <- Vectorize(.pAbeforeB,
                        vectorize.args=c("dateA","dateB", "muA", "muB", "LA", "LB")) ## end .pAbeforeB


##################
## .pAbeforeBfast
##################
##
## faster version, same mu and length for both sequences
## already vectorised for dateA and dateB
.pAbeforeBfast <- function(dateA, dateB, mu, L, maxNbDays=100){
    ## proba dist for both haplo
    temp <- .dTimeSeq(mu, L, maxNbDays)
    days <- temp[[1]]
    p <- temp[[2]]/sum(temp[[2]]) # scale to proba mass function

    ## days for A and B
    nbDays <- as.integer(round(difftime(dateB,dateA,units="days"))) # dateA - dateB, in days

    ## function for one comparison
    f1 <- function(Dt,max){ # does not work for Dt < 0 (have to reverse proba after)
        if(is.na(Dt)) return(NA)
        if(Dt>max) return(1)
        if(round(Dt)==0){
            temp <- sapply(1:(max-1), function(i) p[i]*sum(p[(i+1):max]))
            return(sum(temp))
        }
        term1 <- sum(p[1:Dt])
        idx <- seq(2,by=1,length=(max-Dt))
        temp <- sapply(idx, function(i) sum(p[i:max]))
        term2 <- sum( p[(Dt+1):max] * temp)

        return(term1+term2)
    }

    ## computations
    distribSize <- length(days)
    res <- sapply(nbDays, f1, max=distribSize)
    res[nbDays<0] <- 1-res[nbDays<0] # reverse proba for negative time diff

    return(res)
} # end .pAbeforeBfast







#######################################################
#######################################################

########################
## seqTrack - basic version
########################
##
## - x is a matrix giving weights x[i,j] such that:
## 'i is ancestor j'
## - prox.mat is a directed proximity measure, so that prox.mat[i,j] is
## the 'proximity when going from i to j'
##
seqTrack.default <- function(x, x.names, x.dates, best=c("min","max"),
                     prox.mat=NULL, mu=NULL, haplo.length=NULL, ...){

    ## CHECKS ##
    best <- match.arg(best)
    if(best=="min"){
        best <- min
        which.best <- which.min
    } else {
        best <- max
        which.best <- which.max
    }

    if(length(x.names) != length(x.dates)){
        stop("inconsistent length for x.dates")
    }

    if(is.character(x.dates)){
        msg <- paste("x.dates is a character vector; " ,
                     "please convert it as dates using 'as.POSIXct'" ,
                     "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
        stop(msg)
    }

    x <- as.matrix(x)

    if(!is.null(prox.mat) && !identical(dim(prox.mat),dim(x))) {
        stop("prox.mat is provided but its dimensions are inconsistent with that of x")
    }

    N <- length(x.names)
    id <- 1:N

    x.dates <- as.POSIXct(round.POSIXt(x.dates,units="days")) # round dates to the day

    temp <- as.vector(unique(x))
    D.ARE.MUT <- all(temp-round(temp,10)<1e-14)


    ## rename dimensions using id
    colnames(x) <- rownames(x) <- id
    if(!is.null(prox.mat)){
        colnames(prox.mat) <- rownames(prox.mat) <- id
    }

    if(length(x.names) != nrow(x)){
        stop("inconsistent dimension for x")
    }



    ## AUXILIARY FUNCTIONS ##
    ## test equality in floats
    test.equal <- function(val,vec){
        return(abs(val-vec) < 1e-12)
    }


    ## return the names of optimal value(s) in a named vector
    which.is.best <- function(vec){
        res <- names(vec)[test.equal(best(vec), vec)]
        return(res)
    }


    ## select among different possible ancestors
    selAmongAncestors <- function(idx,ances){
        ## Choose the most connected ancestor, given prox.mat
        if(!is.null(prox.mat)){ # if we've got no other info
            toKeep <- test.equal(max(prox.mat[ances,idx]), prox.mat[ances,idx])
            ances <- ances[toKeep]
        }

        ## If several ancestors remain, take the one closest to the average generation time.
        if(length(ances)>1){
            if(!D.ARE.MUT | is.null(mu) | is.null(haplo.length)) { # if we don't have mutation rates / haplo length, or if dist. are not nb of mutations
                ances <- ances[which.min(x.dates[ances])] # take the oldest ancestor
            } else { # if distances are mutations and we've got mu and L
                timeDiff <- as.numeric(difftime(x.dates[idx], x.dates[ances], units="day")) # days between candidates and target
                ##nbGen <- round(timeDiff / gen.time) # number of generations
                nbMut <- x[ances, idx]
                prob <- dbinom(nbMut, timeDiff*haplo.length, mu)
                ances <- ances[which.max(prob)] # take the most likely ancestor
            }
        }

        return(ances)
    }


    ## findAncestor
    findAncestor <- function(idx){ # returns the index of one seq's ancestor
        candid <- which(x.dates < x.dates[idx])
        if(length(candid)==0) return(list(ances=NA, weight=NA))
        if(length(candid)==1) return(list(ances=candid, weight=x[candid, idx]))
        ancesId <- as.numeric(which.is.best(x[candid, idx]))
        if(length(ancesId)>1) {
            ancesId <- selAmongAncestors(idx,ancesId) # handle several 'best' ancestors
        }
        return(list(ances=ancesId, weight=x[ancesId, idx])) # Id of the ancestor
    }


    ## BUILD THE OUTPUT ##
    res <- sapply(id, findAncestor)
    res <- data.frame(ances=unlist(res[1,]), weight=unlist(res[2,]))
    ances.date <- x.dates[res[,1]]
    res <- cbind.data.frame(id,res, date=x.dates, ances.date)
    rownames(res) <- x.names

    class(res) <- c("seqTrack","data.frame")

    return(res)
} # end seqTrack.default

#######################################################
#######################################################










## ##############
## ## seqTrackG - graph version of SeqTrack
## ##############
## ##
## ## - x is a matrix giving weights x[i,j] such that:
## ## 'i is ancestor j'; the algo looks for maximal weight branching
## ##
## ## - prox.mat is a directed proximity measure, so that prox.mat[i,j] is
## ## the 'proximity when going from i to j'
## ##
## seqTrackG.default <- function(x, x.names, x.dates, best=c("min","max"), force.temporal.order=TRUE,
##                               res.type=c("seqTrack", "graphNEL"), ...){

##     ## CHECKS ##
##     if(!require("graph")) stop("the graph package is not installed")
##     if(!require("RBGL")) stop("the RBGL package is not installed")
##     if(!exists("edmondsOptimumBranching")) {
##         stop("edmondsOptimumBranching does not exist; \nmake sure to use the latest Bioconductor (not CRAN) version of RBGL")
##         cat("\nWould you like to try and install latest version of RBGL (needs internet connection)\n y/n: ")
##         ans <- tolower(as.character(readLines(n = 1)))
##         if(ans=="y"){
##             source("http://bioconductor.org/biocLite.R")
##             biocLite("RBGL")
##         }
##     }

##     best <- match.arg(best)
##     res.type  <- match.arg(res.type)

##     if(length(x.names) != length(x.dates)){
##         stop("inconsistent length for x.dates")
##     }

##     if(is.character(x.dates)){
##         msg <- paste("x.dates is a character vector; " ,
##                      "please convert it as dates using 'as.POSIXct'" ,
##                      "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
##         stop(msg)
##     }

##     x <- as.matrix(x)

##     x.dates <- as.POSIXct(round.POSIXt(x.dates,units="days")) # round dates to the day

##     if(length(x.names) != nrow(x)){
##         stop("inconsistent dimension for x")
##     }


##     ## HANDLE BEST==MIN ##
##     ## reverse x by a translation
##     x.old <- x
##     areZero <- x<1e-14
##     x <- (max(x)+1) - x
##     x[areZero] <- 0
##     diag(x) <- 0


##     ## BUILD THE GRAPH ##
##     ## make prox matrix with temporal consistency
##     if(force.temporal.order){
##         x[outer(myDates, myDates, ">=")] <- 0
##     }

##     ## tweak to get around a bug in as(...,"graphNEL") - looses edge weights
##     ## to replace with commented line once fixed in CRAN release
##     ## myGraph <- as(x, "graphNEL")
##     myGraph <- as(new("graphAM", x, edgemode="directed", values=list(weight=1)), "graphNEL")


##     ## CALL EDMONDSOPTIMUMBRANCHING  ##
##     temp <- edmondsOptimumBranching(myGraph)


##     ## SHAPE OUTPUT ##
##     if(res.type=="seqTrack"){
##         ## reorder output from edmondsOptimumBranching
##         N <- length(x.names)
##         myLev <- x.names
##         ances <- as.integer(factor(temp$edgeList["from",], levels=myLev))
##         desc <- as.integer(factor(temp$edgeList["to",], levels=myLev))
##         newOrd <- order(desc)
##         desc <- 1:N
##         ances <- ances[newOrd]
##         weights <- as.numeric(temp$weights)[newOrd]

##         ## create the data.frame
##         res <- data.frame(id=1:N)
##         hasNoAnces <- difftime(myDates, min(myDates), units="secs") < 1 # 1 sec resolution for dates
##         res$weight <- res$ances <- 1:N
##         res$date <- x.dates

##         ## fill in the d.f. with correct values
##         res$ances[!hasNoAnces] <- ances
##         res$ances[hasNoAnces] <- NA

##         res$weight[!hasNoAnces] <- weights
##         res$weight[hasNoAnces] <- NA

##         res$ances.date <- x.dates[res$ances]

##         res[is.na(res)] <- NA # have clean NAs
##         row.names(res) <- x.names
##         class(res) <- c("seqTrack","data.frame")

##         ## handle best==min
##         if(best=="min"){
##             res$weight <- max(x.old)+1 - res$weight
##         }


##     }

##     if(res.type=="graphNEL"){
##         ## handle optim==min
##         if(best=="min"){
##             temp$weights <- max(x.old)+1 - temp$weights
##         }
##         res <- ftM2graphNEL(t(temp$edgeList), W=temp$weights, edgemode="directed")
##     }


##     ## RETURN RESULT ##
##     return(res)

## } # end seqTrackG







################
## plotSeqTrack
################
plotSeqTrack <- function(x, xy, useArrows=TRUE, annot=TRUE, labels=NULL, dateRange=NULL,
                         col=NULL, bg="grey", add=FALSE, quiet=FALSE,
                         support=NULL, support.thres=0.5, timeorder.thres=NULL,
                         mu=NULL, seq.length=NULL,
                         col.pal=heat.colors, plot=TRUE,...){

    ## CHECKS ##
    if(!inherits(x,"seqTrack")) stop("x is not a seqTrack object")
    ##if(ncol(x) != 5) stop("x does not have five columns")
    if(ncol(xy) != 2) stop("xy does not have two columns")
    if(nrow(xy) != nrow(x)) stop("x and xy have inconsistent dimensions")
    if(!is.null(timeorder.thres) & (is.null(mu) | is.null(seq.length)) ){
        stop("timeorder.thres provided without mu and seq.length.")
    }
    if(!is.null(support)){
        if(length(support)!=nrow(xy)) stop("Inconsistent length for support.")
    }

    isAmbig <- NULL

    ## RECYCLE COL
    if(!is.null(col)){
        col <- rep(col,length=nrow(x))
    }

    ## HANDLE COLOR PALETTE
    if(is.null(col) & !is.null(support)){ # use palette iff support provided without col
        opal <- palette()
        on.exit(palette(opal))
        palette(col.pal(100))
    }

    ## DEFAULT LABELS
    if(is.null(labels)){
        labels <- rownames(x)
    }

    ## SUBSET DATA (REMOVE NAs) ##
    isNA <- is.na(x[,2])
    x <- x[!isNA,,drop=FALSE]
    xy.all <- xy ## used to retrieve all coordinates
    xy <- xy[!isNA,,drop=FALSE]
    if(!is.null(labels)){ # subset labels
        labels <- labels[!isNA]
    }
    if(!is.null(col)){ # subset colors
        col <- col[!isNA]
    }
    if(!is.null(support)){
        support <- support[!isNA] # subset support
    }


    ## FIND AMBIGUOUS ANCESTRIES ##
    if(!is.null(support)){
        isAmbig <- support < support.thres
    }


    ## FIND AMBIGUOUS TEMPORAL ORDERING ##
    if(!is.null(timeorder.thres)){
        temp <- .pAbeforeBfast(x$ances.date, x$date, mu, seq.length)
        if(is.null(isAmbig)){
            isAmbig <- temp < timeorder.thres
        } else {
            isAmbig <- isAmbig | (temp < timeorder.thres)
        }
    }


    ## FIND SEGMENTS COORDS ##
    from <- unlist(x[,2])
    to <- unlist(x[,1])

    x.from <- xy.all[from,1]
    y.from <- xy.all[from,2]
    x.to <- xy.all[to,1]
    y.to <- xy.all[to,2]


    ## FIND THE COLOR FOR EDGES ##
    if(is.null(col)){
        if(!is.null(support)){
            opalette <- palette()
            on.exit(palette(opalette))

            w <- support/max(support,na.rm=TRUE) # support from 0 to 1
            w <- 1 + ((1-w)*99)
            col <- w
        } else{
            col <- rep("black", length(from))
        }
    }


    ## THIS WAS USED WHEN COLOR REPRESENTED THE NUMBER OF MUTATIONS ##
    ##  if(is.null(col)){
    ##         w <- as.numeric(x[,3])
    ##         w <- max(w) - w
    ##         w <- w-min(w)
    ##         w <- 1+ w/max(w) * 99

    ##         opalette <- palette()
    ##         on.exit(palette(opalette))
    ##         palette(heat.colors(100))

    ##         col <- w
    ##     }


    ## HANDLE RANGE OF DATES ##
    if(!is.null(dateRange)){

        if(is.character(dateRange)){
            msg <- paste("dateRange is a character vector; " ,
                         "please convert it as dates using 'as.POSIXct'" ,
                         "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
            stop(msg)
        }

        if(any(is.na(dateRange))){
            stop("NA in dateRange")
        }

        dates <- x$date
        toKeep <- (dates > min(dateRange)) & (dates < max(dateRange))
        if(sum(toKeep)==0) {
            if(!quiet) cat("\nNo item in the specified date range.\n")
            return(NULL)
        }

        ## do the subsetting
        x.from <- x.from[toKeep]
        y.from <- y.from[toKeep]
        x.to <- x.to[toKeep]
        y.to <- y.to[toKeep]
        col <- col[toKeep]
        xy <- xy[toKeep,,drop=FALSE]
        x <- x[toKeep,,drop=FALSE]
        if(!is.null(isAmbig)) { # subset isAmbig
            isAmbig <- isAmbig[toKeep]
        }
        if(!is.null(labels)){ # subset labels
            labels <- labels[toKeep]
        }
    }



    ## DO THE PLOTTING ##
    if(plot){
        obg <- par("bg")
        on.exit(par(bg=obg))
        if(!add){
            par(bg=bg)
            plot(xy, type="n")
        }
    }

    ## ARROWS
    if(useArrows & plot){
        ## handle segments/arrows with length 0 ##
        nullLength <- (abs(x.from-x.to)<1e-10) & (abs(y.from-y.to)<1e-10)

        if(any(isAmbig)){ # plot arrows & segments
            arrows(x.from[!isAmbig & !nullLength], y.from[!isAmbig & !nullLength],
                                    x.to[!isAmbig & !nullLength], y.to[!isAmbig & !nullLength],
                                    col=col[!isAmbig & !nullLength], angle=15, ...)
            segments(x.from[isAmbig], y.from[isAmbig],
                     x.to[isAmbig], y.to[isAmbig], col=col[isAmbig],...)
        } else{ # plot all arrows
            arrows(x.from[!nullLength], y.from[!nullLength],
                                    x.to[!nullLength], y.to[!nullLength],
                   col=col[!nullLength], angle=15, ...)
        }
    } else{
        ## SEGMENTS
        if(plot) segments(x.from, y.from, x.to, y.to, col=col,...)
    }


    if(annot & plot) {
        text(xy,lab=labels, font=2)
    }

    if(any(nullLength) & plot) {
        sunflowerplot(x.from[nullLength], y.from[nullLength], seg.lwd=2, size=1/6,
                      col=col[nullLength], seg.col=col[nullLength], add=TRUE, ...)
        points(x.from[nullLength], y.from[nullLength], col=col[nullLength], cex=2, pch=20, ...)
    }


    ## RESULT ##
    res <- data.frame(x.from, y.from, x.to, y.to, col=col)
    if(!is.null(isAmbig)) {
        res <- cbind.data.frame(res, isAmbig)
    }
    return(invisible(res))
} # end plotSeqTrack












#####################
## optimize.seqTrack
#####################
##
## TODO:
## 1) Change the output to retain xxx simulations | ok. -- done.
## 2) VECTORIZE mu and seq.length, recycle if needed with a warning
## 3) uncomment, adapt, and test code for missing data
##
## optimize.seqTrack.default <- function(x, x.names, x.dates, typed.chr=NULL, mu=NULL, seq.length=NULL,
##                                       thres=0.2, best=c("min","max"), prox.mat=NULL, nstep=10, step.size=1e3,
##                                       rDate=.rTimeSeq, arg.rDate=NULL, rMissDate=.rUnifDate, ...){


##     ## CHECKS ##
##     best <- match.arg(best)
##     if(best=="min"){
##         which.best <- which.min
##     } else {
##         which.best <- which.max
##     }

##     if(length(x.names) != length(x.dates)){
##         stop("inconsistent length for x.dates")
##     }

##     if(is.character(x.dates)){
##         msg <- paste("x.dates is a character vector; " ,
##                      "please convert it as dates using 'as.POSIXct'" ,
##                      "\n(making sure dates are given as 'YYYY/MM/DD' or 'YYYY-MM-DD').", sep="")
##         stop(msg)
##     }

##     isMissDate <- is.na(x.dates)

##     if(!identical(rDate, .rTimeSeq)){
##         if(is.null(arg.rDate)){
##             warning("Specific time distribution specified without arguments.")
##             arg.rDate <- list(n=step.size)
##         } else {
##             if(!is.list(arg.rDate)) stop("If provided, arg.rDate must be a list.")
##             if(!is.null(arg.rDate$n)) {
##                 warning("arg.rDate$n is provided, but will be replaced by step.size.")
##             }
##             arg.rDate$n <- step.size
##         }
##     }


##     N <- length(x.names)
##     id <- 1:N
##     ## if(length(mu) < N) { # recycle mu
##     ##         mu <- rep(mu, length=N)
##     ##     }
##     ##     if(length(seq.length) < N) {# recycle seq.length
##     ##         seq.length <- rep(seq.length, length=N)
##     ##     }


##     ## handle typed.chr, mu, seq.length
##     if(identical(rDate, .rTimeSeq)){
##         if(is.null(typed.chr)|is.null(mu)|is.null(seq.length)){
##             stop("typed.chr, mu, and seq.length must be provided if rDate is .rTimeSeq")
##         }

##         if(!is.list(typed.chr)) {
##             stop("typed.chr must be a list")
##         }
##         if(length(typed.chr)!=N) {
##             stop("typed.chr has an inconsistent length")
##         }

##         if(is.null(names(mu))) stop("mu has no names")
##         if(is.null(names(seq.length))) stop("seq.length has no names")
##         if(any(mu > 1)) stop("mu has values > 1")
##         if(any(mu < 0)) stop("mu has negative values")

##         if(!identical(names(mu) , names(seq.length))) stop("Names of mu and seq.length differ.")
##         if(any(!unique(unlist(typed.chr)) %in% names(mu))) {
##             stop("Some chromosomes indicated in typed.chr are not in mu.")
##         }

##         list.mu <- lapply(typed.chr, function(e) mu[e])
##         list.seq.length <- lapply(typed.chr, function(e) seq.length[e])
##     }

##     x.dates <- as.POSIXct(round.POSIXt(x.dates,units="days")) # round dates to the day

##     x <- as.matrix(x)

##     if(!is.null(prox.mat) && !identical(dim(prox.mat),dim(x))) {
##         stop("prox.mat is provided but its dimensions are inconsistent with that of x")
##     }


##     ## rename dimensions using id
##     colnames(x) <- rownames(x) <- id

##     if(length(x.names) != nrow(x)){
##         stop("inconsistent dimension for x")
##     }


##     ## SET THRESHOLD IF NEEDED ## ## NO LONGER USED
##     ##  if(is.null(thres)){
##     ##         thres <- sum(seqTrack(x.names=x.names, x.dates=x.dates, W=W,
##     ##                                     best=best, prox.mat=prox.mat, ...)$weight, na.rm=TRUE)
##     ##     }


##     ## AUXILIARY FUNCTIONS ##

##     ## to compare results -> returns a list of length two: logical, and the value of the res
##     val.res <- function(res){
##         return(sum(res$weight, na.rm=TRUE))
##     }


##     ## DO THE OPTIMISATION ##
##     RANGE.DATES <- as.integer(round(diff(range(x.dates, na.rm=TRUE)))) # time window of the sample, in days
##     NB.DATES.TO.SIM <- sum(!isMissDate)


##     ## for loop is not to slow for < 1e6 rep
##     ## and allows not to handle huge objects
##     ## (which would grow exponentially)

##     ##    res.best <- res.ini # initialization


##     ## DEFINE OUTPUTS ##
##     ances <- integer(0)
##     date <- character(0)
##     ances.date <- character(0)
##     valRes <- numeric(0)


##     ## DEFAULT CASE: NO MISSING DATES
##     if(!any(isMissDate)){
##         ## dates initialisation, taken from initial prior
##         ## If dates distrib is .rTimeSeq
##         if(identical(rDate, .rTimeSeq)){
##             newDates <- sapply(1:N, function(i)
##                                rDate(n=step.size, mu=list.mu[[i]], L=list.seq.length[[i]],
##                                      maxNbDays=RANGE.DATES))
##         } else { ## Else, any other distrib with free arguements
##             newDates <- sapply(1:N, function(i) do.call(rDate, arg.rDate))
##         }

##         newDates <- t(newDates)*24*3600 + x.dates

##         ## >> one step of 'step.size' simulations, all with same prior << ##
##         for(i in 1:nstep){
##             ## >> each step contains 'step.size' iterations << ##
##             for(j in 1:step.size){
##                 myDates <- as.POSIXct(newDates[,j])

##                 res.new <- seqTrack(x, x.names=x.names, x.dates=myDates,
##                                     best=best, prox.mat=prox.mat, ...)

##                 ##ances <- cbind(ances, res.new$ances) # not needed now
##                 date <- cbind(date, as.character(res.new$date))
##                 ##ances.date <- cbind(ances.date, as.character(res.new$ances.date)) # not needed now
##                 valRes <- c(valRes, val.res(res.new))
##                 ##}
##             } # end for j

##             ## retain a given % (thres) of the dates ##
##             toKeep <- valRes <= quantile(valRes, thres) ## NOT WORKING FOR optim==max !!!
##             valRes <- valRes[toKeep]

##             date <- date[,toKeep,drop=FALSE] # retained posterior

##             ## DEBUGING ##
##             ## cat("\ntoKeep:\n")
##             ##             print(toKeep)
##             ##             cat("\nhead date (posterior):\n")
##             ##             print(head(date))
##             ## END DEBUGING ##

##             newDates <- apply(date, 1, function(vec)
##                               sample(vec, size=step.size, replace=TRUE)) # new prior
##             newDates <- t(newDates)

##             ## stop if all dates are fixed
##             if(all(apply(newDates, 1, function(r) length(unique(r))==1))){
##                 cat("\nConvergence reached at step",i,"\n")
##                 break # stop the algorithm
##             }

##             ## re-initialize posterior distributions
##             if(i<nstep){
##                 ## ances <- integer(0) # not needed now
##                 date <- character(0)
##                 ## ances.date <- character(0) # not needed now
##                 valRes <- numeric(0)
##             } # end if

##         } # end for i

##         ##  ## dates: new prior taken from obtained posterior
##         ##             if(length(valRes)==0) { # if no simul are retained
##         ##                 warning(paste("No simulation was retained at the given threshold at step",i))
##         ##             } else {
##         ##  if(optim=="min"){ # define weights for further samplings
##         ##                     w <- max(valRes,na.rm=TRUE) - valRes
##         ##                     w <- w/sum(w)
##         ##                 } else {
##         ##                     w <- valRes
##         ##                     w <- w/sum(w)
##         ##                 }

##         ## newDates <- apply(date, 1, function(vec) #  used a weighted sampling
##         ##                                  sample(vec, size=step.size, replace=TRUE, prob=w))


##     } # end if(!any(isMissDate))


##     ##  ## OTHER CASE: HANDLE MISSING DATES
##     ##     if(any(isMissDate)){

##     ##         ## Handle distribution and its parameters ##
##     ##         argList <- list(...)

##     ##         if(is.null(argList$dateMin) & identical(rMissDate, .rUnifDate)){ # earliest date
##     ##             argList$dateMin <- min(x.dates,na.rm=TRUE)
##     ##         } else {
##     ##             argList$dateMin[is.na(argList$dateMin)] <- min(x.dates,na.rm=TRUE)
##     ##         }
##     ##         if(is.null(argList$dateMax) & identical(rMissDate, .rUnifDate)){ # latest date
##     ##             argList$dateMax <- max(x.dates,na.rm=TRUE)
##     ##         } else {
##     ##             argList$dateMax[is.na(argList$dateMax)] <- max(x.dates,na.rm=TRUE)
##     ##         }

##     ##         argList$n <- sum(isMissDate)

##     ##         ## Make simulations ##
##     ##         for(i in 1:nstep){
##     ##             myDates <- x.dates
##     ##             ## distribution for available dates
##     ##             myDates[!isMissDate] <- myDates[!isMissDate] +
##     ##                 rDate(n=NB.DATES.TO.SIM, mu=mu, L=seq.length, maxNbDays=2*RANGE.DATES)
##     ##             ## distribution for missing dates
##     ##             myDates[isMissDate] <- do.call(rMissDate, argList)

##     ##             res.new <- seqTrack(x.names=x.names, x.dates=myDates, W=W, optim=optim, prox.mat=prox.mat, ...)

##     ##             valRes[i] <- sum(res.new$weight,na.rm=TRUE)
##     ##             if(use.new.res(res.best, res.new)){
##     ##                 res.best <- res.new
##     ##             }
##     ##         }
##     ##     }


##     ## RESULT ##

##     ## reconstruct the result with new dates
##     res <- lapply(1:ncol(date), function(i)
##                    seqTrack(x=x, x.names=x.names, x.dates=as.POSIXct(date[,i]),
##                                     best=best, prox.mat=prox.mat, ...))
##     ances <- data.frame(lapply(res, function(e) e$ances))
##     ances <- matrix(as.integer(unlist(ances)), nrow=nrow(ances))

##     ances.date <- data.frame(lapply(res, function(e) as.character(e$ances.date)))
##     ances.date <- matrix(as.character(unlist(ances.date)), nrow=nrow(ances.date))

##     res <- list(ances=ances, date=date, ances.date=ances.date, valsim=valRes)
##     return(res)

## } # optimize.seqTrack






## #################
## ## get.result.by
## #################
## .get.result.by <- function(x, ...){
##     dat <- list(...)
##     if(length(dat)==0) return(x)


##     ## DEFINE NEW VALUES ##

##     convertElem <- function(e){
##         if(class(e)=="DNAbin") {
##             e <- as.matrix(e)
##         }
##         ori.dim <- dim(e)
##         e <- as.character(e)
##         dim(e) <- ori.dim
##         return(e)
##     }


##     dat <- lapply(dat,convertElem)
##     dat <- as.matrix(data.frame(dat))

##     newval <- apply(dat, 1, function(vec) paste(vec, collapse=""))
##     newval <- unclass(factor(newval))
##     newlev <- levels(newval)


##     ## if x is a single output of seqTrack
##     if(is.vector(x$ances)){
##         newId <- newval # new values
##         newAnces <- newval[x$ances] # new values
##         ## make output
##         res <- x
##         res$id <- newId
##         res$ances <- newAnces
##         attr(res$ances, "levels") <- newlev
##     }


##     ## if x is an optimize.seqTrack output
##     if(is.matrix(x$ances)){
##         res <- x
##         ori.ncol <- ncol(res$ances)
##         res$ances <- matrix(newval[res$ances], ncol=ori.ncol)
##         attr(res$ances, "levels") <- newlev
##     }

##     ## method for haploGen
##     if(inherits(x,"haploGen")){
##         res <- x
##         ances.id <- match(x$ances, labels(x))
##     }

##     return(res)

## } # end get.result.by






## #################
## ## get.consensus
## #################
## .get.consensus <- function(orires, listres, mode=c("majority","best")){
##     ## handle mode
##     mode <- match.arg(mode)

##     res <- orires

##     if(mode=="majority"){
##         nbDraws <- 0

##         ## tables of occurences of ancestors
##         temp <- apply(listres$ances, 1, table)

##         ## compute compromise
##         if(!is.list(temp)){
##             newances <- temp
##             ances.support <- rep(1,length(temp))
##         } else {
##             f1 <- function(tab){
##                 if(length(tab)==0) return(NA)

##                 res <- names(tab)[tab==max(tab)]
##                 ## if(length(res)==1) return(res)
##                 ##             return(NA)
##                 if(length(res)>1) {
##                     nbDraws <- nbDraws+1
##                 }
##                 return(res[1])
##             }

##             newances <- sapply(temp, f1)
##             ances.support <- sapply(temp, function(e) max(e, na.rm=TRUE)/sum(e, na.rm=TRUE))
##             ances.support[is.na(newances)] <- NA
##         }

##         ## form the output
##         olev <- levels(orires$ances)
##         res$ances <- newances
##         levels(res$ances) <- olev
##         res$support <- ances.support
##         res$weight <- rep(1, length(res$date))

##         if(is.numeric(listres$ances)){
##             res$ances <- as.numeric(res$ances)
##         }
##         cat("\nThere were\n",nbDraws, "draws.\n")
##     } # end majority


##     if(mode=="best"){
##         toKeep <- which.min(listres$valsim)
##         nbDraws <- sum(listres$valsim < (min(listres$valsim) + 1e-10 )) -1
##         cat("\nThere were\n",nbDraws, "draws.\n")

##         res$ances <- listres$ances[,toKeep]
##         res$inf.date <- listres$date[,toKeep]
##         res$ances.date <- listres$ances.date[,toKeep]
##         res$weight <- rep(-1, length(res$date))
##     }

##     return(res)
## } # end get.consensus







###########################
## get.likelihood.seqTrack
###########################
get.likelihood.seqTrack <-function(x, method=("genetic"), mu=NULL, seq.length=NULL,...){
    method <- match.arg(method)
    if(method=="genetic"){ # p(nb mutations occur in the time interval)
        if(any(na.omit(x$weight - round(x$weight)) > 1e-10)){
            warning("Non-integer weights: number of mutations expected in x$weight.")
        }

        if(is.null(mu)) stop("mu is required.")
        if(is.null(seq.length)) stop("seq.length is required.")

        dates <- as.POSIXct(x$date)
        anc.dates <- as.POSIXct(x$ances.date)
        nb.days <- abs(as.integer(anc.dates-dates))
        nb.mut <- x$weight
        ##mu <- mu/365
        ##mu <- mu*nb.days

        res <- dbinom(nb.mut, size=seq.length*nb.days, prob=mu)
    } else{
        cat("Method not implemented.")
    }

    return(res)
} # end get.likelihood.seqTrack






###############
## seqTrack.ml
###############
## seqTrack.ml <- function(aln, x.dates, best=c("min","max"), ...){

## } # end seqTrack.ml









##########################
## as("seqTrack", "graphNEL")
##########################
setOldClass("seqTrack")
setAs("seqTrack", "graphNEL", def=function(from){
    if(!require(ape)) stop("package ape is required")
    if(!require(graph)) stop("package graph is required")

    ori.labels <- rownames(from)
    from <- from[!is.na(from$ances),,drop=FALSE]


     ## CONVERT TO GRAPH
    res <- ftM2graphNEL(ft=cbind(ori.labels[from$ances], ori.labels[from$id]), W=from$weight, edgemode = "directed", V=ori.labels)
    return(res)
})

