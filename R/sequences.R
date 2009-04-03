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
    x <- x[,toKeep]

    ## build output
    res <- df2genind(x, pop=pop, ploidy=1, ncode=1, type="codom")
    res$call <- match.call()

    return(res)
} # end DNAbin2genind
