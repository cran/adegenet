######################################
##
## The code below implements import
## from alignement data.
##
######################################



################
# DNAbin2genind
################
DNAbin2genind <- function(x, pop=NULL, exp.char=c("a","t","g","c"), na.char=NULL, polyThres=1/100){

    ## misc checks
    if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
    if(!require(ape)) stop("The package ape is required.")

    ## DNA bin to matrix of characters
    x <- as.character(x) # should output a matrix

    if(is.list(x)) { # if this is a list
        temp <- unique(sapply(x,length)) # check lengths of sequences
        if(length(temp)>1) stop("Sequences have different length - please use alignements only.")
        else{ # if sequences have same length, build the matrix
            temp <- names(x)
            x <- t(as.data.frame(x))
            rownames(x) <- temp
        }
    }

    if(is.null(colnames(x))) {
        colnames(x) <- 1:ncol(x)
    }

    ## replace NAs
    if(is.null(na.char)){
        if(is.null(exp.char)) stop("both exp.char and na.char are NULL")
        temp <- paste(exp.char, collapse="", sep="")
        if(any(exp.char=="-")) {
            temp <- paste("-",temp, sep="") # string '-' must begin the regexp
        }
        temp <- paste("[^", temp, "]", sep="") # anything but the expected is NA
        x <- gsub(temp,NA,x)
    } else {
        temp <- paste(na.char, collapse="", sep="")
        if(any(na.char=="-")) {
            temp <- paste("-",temp, sep="") # string '-' must start the regexp
        }
        temp <- paste("[", temp, "]", sep="")
        x <- gsub(temp,NA,x)
    }

    ## keep only columns with polymorphism (i.e., SNPs)
    isPoly <- function(vec){
        N <- sum(!is.na(vec)) # N: number of sequences
        temp <- table(vec)/N
        if(sum(temp > polyThres) >= 2) return(TRUE)
        return(FALSE)
    }

    toKeep <- apply(x, 2, isPoly)
    if(sum(toKeep)==0) stop("No polymorphic site detected")
    x <- x[,toKeep]

    ## build output
    res <- df2genind(x, pop=pop, ploidy=1, ncode=1, type="codom")
    res$call <- match.call()

    return(res)
} # end DNAbin2genind







################
# alignment2genind
################
alignment2genind <- function(x, pop=NULL, exp.char=c("a","t","g","c"), na.char="-", polyThres=1/100){

    ## misc checks
    if(!require(seqinr)) stop("The package seqinr is required.")
    if(!inherits(x,"alignment")) stop("x is not a alignment object")
    N <- length(x$seq)
    if(!is.null(x$nam) && length(x$nam)!=N) stop("Inconsistent names in x (length of x$nam and x$seq do not match). ")


    ## convert alignment to matrix of characters
    mat <- sapply(x$seq, s2c, USE.NAMES=FALSE)
    if(nrow(mat)!=x$nb){
        mat <- t(mat)
    }

    rownames(mat) <- x$nam

    if(is.null(colnames(x))) {
        colnames(mat) <- 1:ncol(mat)
    }

    ## replace NAs
    if(is.null(na.char)){
        if(is.null(exp.char)) stop("both exp.char and na.char are NULL")
        temp <- paste(exp.char, collapse="", sep="")
        if(any(exp.char=="-")) {
            temp <- paste("-",temp, sep="") # string '-' must begin the regexp
        }
        temp <- paste("[^", temp, "]", sep="") # anything but the expected is NA
        mat <- gsub(temp,NA,mat)
    } else {
        temp <- paste(na.char, collapse="", sep="")
        if(any(na.char=="-")) {
            temp <- paste("-",temp, sep="") # string '-' must start the regexp
        }
        temp <- paste("[", temp, "]", sep="")
        mat <- gsub(temp,NA,mat)
    }

    ## keep only columns with polymorphism (i.e., SNPs)
    isPoly <- function(vec){
        N <- sum(!is.na(vec)) # N: number of sequences
        temp <- table(vec)/N
        if(sum(temp > polyThres) >= 2) return(TRUE)
        return(FALSE)
    }

    toKeep <- apply(mat, 2, isPoly)
    if(sum(toKeep)==0) stop("No polymorphic site detected")

    mat <- mat[,toKeep, drop=FALSE]

    ## build output
    res <- df2genind(mat, pop=pop, ploidy=1, ncode=1, type="codom")
    res$call <- match.call()

    if(!is.null(x$com)){
        res@other$com <- x$com
    }

    return(res)
} # end alignment2genind




























## ###############
## ## transiProb
## ###############
## ##
## ## proba/distance based on transition prob from one sequence to another
## ## time is taken into account
## ## output: matrix with term proba(rowIdx to colIdx)
## ##
## transiProb <- function(x, mu, dates, result=c("prob","dist")){
##     ## MISC CHECKS ##
##     if(!inherits(x,"DNAbin")) stop("x is not a DNAbin object")
##     if(!require(ape)) stop("The package ape is required.")
##     result <- match.arg(result)

##     ## COMPUTATIONS ##

##     ## get numbers of differing nucleotides between sequences
##     seq.length <- ncol(as.matrix(x))
##     D <- as.matrix(dist.dna(x, model="raw")) * seq.length
##     ## if(sum(D-round(D)) > 1e-10){ # make sure we've got integers there
##     ##         warning("Number of nucleotides are not all integers")
##     ##     }
##     D <- round(D)

##     ## compute matrix T (time between sequences)
##     if(inherits(dates,"POSIXct")){ # dates in POSIXct format
##         temp <- outer(dates, dates, difftime, unit="days")
##         T <- -matrix(as.numeric(temp),ncol=length(dates))
##     } else { # dates are numeric
##         T <- -outer(dates, dates, "-")
##     }

##     ## spot negative times
##     toSetToNull <- T < 1e-15

##     ## compute proba(no change @ a site) term
##     mu <- mu/365 # express mu per day
##     p1 <- exp(-T*mu) + (1-exp(-T*mu))/4
##     p1[toSetToNull] <- 0
##     res <- dbinom(D, size=seq.length, prob=(1-p1))

##     ## PROCESS/RETURN RESULT
##     if(result=="prob"){ # return probabilities
##         res[toSetToNull] <- 0
##         diag(res) <- 1
##     } else { # return d = -log(proba)
##         res <- -log(res)
##         res[toSetToNull] <- 1e15
##         diag(res) <- 0
##     }

##     return(res)
## } # end transiProb
