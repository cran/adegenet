###################################################################
## Fonctions designed to import files from other softwares
## into genind objects
##
## currently supported formats are :
## .gtx (GENETIX)
## .dat (Fstat)
## .gen (Genepop)
## .stru (STRUCTURE)
##
## Thibaut Jombart, avril 2006
## t.jombart@imperial.ac.uk
##
##################################################################



#####################
# Function df2genind
#####################
df2genind <- function(X, sep=NULL, ncode=NULL, ind.names=NULL, loc.names=NULL, pop=NULL, missing=NA, ploidy=2, type=c("codom","PA")){

    if(is.data.frame(X)) X <- as.matrix(X)
    if (!inherits(X, "matrix")) stop ("X is not a matrix")

    res <- list()
    type <- match.arg(type)
    ## checkType(type)


    ## type-independent stuff ##
    n <- nrow(X)
    nloc <- ncol(X)
    ploidy <- as.integer(ploidy)
    if(ploidy < 1L) stop("ploidy cannot be less than 1")

    if(is.null(ind.names)) {ind.names <- rownames(X)}
    if(is.null(loc.names)) {loc.names <- colnames(X)}

    ## pop optionnelle
    if(!is.null(pop)){
        if(length(pop)!= n) stop("length of factor pop differs from nrow(X)")
        pop <- as.factor(pop)
    }


    ## PA case ##
    if(toupper(type)=="PA"){
        ## preliminary stuff
        mode(X) <- "numeric"
        rownames(X) <- ind.names
        colnames(X) <- loc.names

        ## Erase entierely non-typed loci
        temp <- apply(X,2,function(c) all(is.na(c)))
        if(any(temp)){
            X <- X[,!temp]
            warning("entirely non-type marker(s) deleted")
        }

        ## Erase entierely non-type individuals
        temp <- apply(X,1,function(r) all(is.na(r)))
        if(any(temp)){
            X <- X[!temp,,drop=FALSE]
            pop <- pop[!temp]
            warning("entirely non-type individual(s) deleted")
        }

        ## erase non-polymorphic loci
        temp <- apply(X, 2, function(loc) length(unique(loc[!is.na(loc)]))==1)
        if(any(temp)){
            X <- X[,!temp,drop=FALSE]
            warning("non-polymorphic marker(s) deleted")
        }

        prevcall <- match.call()

        res <- genind( tab=X, pop=pop, prevcall=prevcall, ploidy=ploidy, type="PA")

        return(res)
    } # end type PA


    ## codom case ##

    ## make sure X is in character mode
    mode(X) <- "character"


    ## find or check the number of coding characters, 'ncode'
    if(is.null(sep)){
        if(!is.null(ncode)) {
            temp <- nchar(X[!is.na(X)])
            if(ncode <  max(temp) ) stop("some character strings exceed the provided ncode.")
        }
        if(is.null(ncode)) {
            temp <- nchar(X[!is.na(X)])
            ncode <- max(temp)
        }
        if((ncode %% ploidy)>0) stop(paste(ploidy,"alleles cannot be coded by a total of",
                                           ncode,"characters", sep=" "))
    }

    ## ERASE ENTIRELY NON-TYPE LOCI AND INDIVIDUALS
    tempX <- X
    if(!is.null(sep)) tempX <- gsub(sep,"",X)
    ## turn NANANA, 00000, ... into NA
    tempX <- gsub("^0*$",NA,tempX)
    tempX <- gsub("(NA)+",NA,tempX)

    ## Erase entierely non-typed loci
    temp <- apply(tempX,2,function(c) all(is.na(c)))
    if(any(temp)){
        X <- X[,!temp,drop=FALSE]
        tempX <- tempX[,!temp,drop=FALSE]
        loc.names <- loc.names[!temp]
        nloc <- ncol(X)
        warning("entirely non-type marker(s) deleted")
    }

    ## Erase entierely non-type individuals
    temp <- apply(tempX,1,function(r) all(is.na(r)))
    if(any(temp)){
        X <- X[!temp,,drop=FALSE]
        tempX <- tempX[!temp,,drop=FALSE]
        pop <- pop[!temp]
        ind.names <- ind.names[!temp]
        n <- nrow(X)
        warning("entirely non-type individual(s) deleted")
    }

    n <- nrow(X)
    ## SET NAs IN X
    X[is.na(tempX)] <- NA

    # ind.names <- rownames(X) this erases the real labels
    # note: if X is kept as a matrix, duplicate row names are no problem


    ## function to fill a matrix of char 'M' with the required
    ## number of zero, targetN being the total number of char required
    fillWithZero <- function(M, targetN){
        naIdx <- is.na(M)
        keepCheck <- any(nchar(M) < targetN)
        while(keepCheck){
            mat0 <- matrix("", ncol=ncol(M), nrow=nrow(M))
            mat0[nchar(M) < targetN] <- "0"
            M <-  matrix(paste(mat0, M, sep=""), nrow=nrow(mat0))
            keepCheck <- any(nchar(M) < targetN)
        }

        ## restore NA (otherwise we're left with "NA")
        M[naIdx] <- NA
        return(M)
    }

    ## CHECK STRING LENGTH IF NO SEPARATOR PROVIDED
    if(is.null(sep) | ploidy==as.integer(1)){
        ##     ## now check all strings and make sure they all have 'ncode' characters
        ##         ## NA are temporarily coded as "00", "000" or "000000" to fit the check
        ##         keepCheck <- any(nchar(X) < ncode)
        ##         missAll <- paste(rep("0",ncode/ploidy),collapse="")
        ##         missTyp <- paste(rep("0",ncode),collapse="")
        ##         X[is.na(X)] <- missTyp

        ##         while(keepCheck){
        ##             mat0 <- matrix("", ncol=ncol(X), nrow=nrow(X))
        ##             mat0[nchar(X) < ncode] <- "0"
        ##             X <-  matrix(paste(mat0, X, sep=""), nrow=nrow(mat0))
        ##             keepCheck <- any(nchar(X) < ncode)
        ##         }

        X <- fillWithZero(X,targetN=ncode)

        ## now split X by allele
        splitX <- list()
        for(i in 1:ploidy){
            splitX[[i]] <- substr(X,1,ncode/ploidy)
            X <- sub(paste("^.{",ncode/ploidy,"}",sep=""),"",X)
        }

    } # END CHECK STRING LENGTH WITHOUT SEP


    ## CHECK STRING LENGTH WITH SEPARATOR PROVIDED
    if(!is.null(sep)){
        if(ploidy > 1){
            temp <- t(as.matrix(as.data.frame(strsplit(X,sep))))
            splitX <- list()
            for(i in 1:ncol(temp)){
                splitX[[i]] <- matrix(temp[,i], nrow=n)
            } # each matrix of splitX contains typing for 1 allele
        } else {
            splitX <- list()
            splitX[[1]] <- X
        }

        ## get the right ncode
        temp <- unlist(splitX)
        temp <- temp[!is.na(temp)]
        ncode <- max(nchar(temp))*ploidy
        splitX <- lapply(splitX, function(Y) fillWithZero(Y,targetN=ncode/ploidy))
    } # END CHECK STRING LENGTH WITH SEP


    ## AT THIS STAGE, splitX IS A LIST OF MATRICES,
    ## EACH GIVING TYPING FOR AN ALLELE

    ## fetch all possible alleles per locus
    loc.all <- list()
    for(i in 1:nloc){
        temp <- unlist(lapply(splitX,function(e) e[,i]))
        loc.all[[i]] <- sort(unique(temp[!is.na(temp)]))
    }

    names(loc.all) <- loc.names
    ## loc.all is a list whose element are vectors of sorted possible alleles at a locus
    temp <- lapply(1:nloc, function(i) matrix(0,nrow=n,ncol=length(loc.all[[i]]),
       dimnames=list(NULL,loc.all[[i]])) )

    names(temp) <- loc.names
    # note: keep rownames as NULL in case of duplicates
    ## temp is a list whose elements are one matrix (indiv x alleles) for each marker

    ## now tables in 'temp' are filled up
    findall <- function(cha,loc.all){
        if(is.na(cha)) return(NULL)
        return(which(cha==loc.all))
    }

    for(k in 1:ploidy){
        for(i in 1:n){
            for(j in 1:nloc){
                allIdx <- findall(splitX[[k]][i,j],loc.all[[j]])
                temp[[j]][i,allIdx] <- temp[[j]][i,allIdx] + 1
                if(is.null(allIdx)) {temp[[j]][i,] <- NA}
            }
        }
    }

    ## beware: colnames are wrong when there is only one allele in a locus
    ## right colnames are first generated
    nall <- unlist(lapply(temp,ncol))
    loc.rep <- rep(names(nall),nall)
    col.lab <- paste(loc.rep,unlist(loc.all,use.names=FALSE),sep=".")

    mat <- as.matrix(cbind.data.frame(temp))
    mat <- mat/ploidy
    colnames(mat) <- col.lab
    rownames(mat) <- ind.names

    if(!is.na(missing)){
      if(missing==0) {mat[is.na(mat)] <- 0}
      if(toupper(missing)=="MEAN") {
        moy <- apply(mat,2,function(c) mean(c,na.rm=TRUE))
        for(j in 1:ncol(mat)) {mat[,j][is.na(mat[,j])] <- moy[j]}
      }
    }

    prevcall <- match.call()

    res <- genind( tab=mat, pop=pop, prevcall=prevcall, ploidy=ploidy, type=type)

    return(res)
} # end df2genind





########################################
# Function read.genetix
# code based on previous ade4 functions
########################################
read.genetix <- function(file=NULL,missing=NA,quiet=FALSE) {
    if(!quiet) cat("\n Converting data from GENETIX to a genind object... \n")


    ## read from file
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed

    if(toupper(.readExt(file)) != "GTX") stop("File extension .gtx expected")
    ## retrieve first infos
    nloc <- as.integer(scan(file,nlines=1,what="character",quiet=TRUE)[1])
    npop <- as.integer(scan(file,nlines=1,skip=1,what="character",quiet=TRUE)[1])
    txt <- scan(file,skip=2,what="character",sep="\n",quiet=TRUE)
    txt <- gsub("\t"," ",txt)
    ## check that nloc is consistent with actual nloc (bug-report 1.2-2.02)
    temp <- temp <- .rmspaces(txt[length(txt)])
    nlocbis <- length(unlist(strsplit(temp, "[[:space:]]+")))-1
    if(nloc != nlocbis) {
        warning(paste("\n== Genetix file error == \n",
                      "Indicated number of locus (", nloc, ")\n",
                      "does not match actual number (", nlocbis, ").\n",
                      "Using ", nlocbis, " as number of locus.\n",
                      "Please check your file.", sep=""))
        nloc <- nlocbis
    }
    loc.names <- txt[seq(1,by=2,length=nloc)]
    txt <- txt[-(1:(nloc*2))]

    ## retrieve populations infos
    pop.names <- vector(mode="character",length=npop)
    pop.nind <- vector(mode="integer",length=npop)
    index <- 1
    temp <- vector(mode="integer",length=npop)
    for(i in 1:npop){
        pop.names[i] <- txt[index]
        pop.nind[i] <- as.numeric(txt[index+1])
        temp[i] <- index
        index <- index + pop.nind[i] + 2
    }
    pop.names <- .rmspaces(pop.names)

    ## retrieve genotypes infos
    txt <- txt[-c(temp,temp+1)]
    txt <- .rmspaces(txt)
    txt <- sapply(1:length(txt),function(i) unlist(strsplit(txt[i],"([[:space:]]+)|([[:blank:]]+)")) )
    X <- t(txt)
    if(ncol(X) == (nloc+1)){
        rownames(X) <- X[,1]
        X <- X[,-1]
    } else{
        rownames(X) <- 1:nrow(X)
    }

    colnames(X) <- loc.names

    ## make a factor "pop" if there is more than one population
    pop <- factor(rep(pop.names,pop.nind))

    ## pass X to df2genind
    res <- df2genind(X=X, ncode=6, pop=pop, missing=missing, ploidy=2)
    res@call <- match.call()

    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end read.genetix



##########################
# Function read.fstat
##########################
read.fstat <- function(file,missing=NA,quiet=FALSE){
    ##if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(toupper(.readExt(file)) != "DAT") stop("File extension .dat expected")

    if(!quiet) cat("\n Converting data from a FSTAT .dat file to a genind object... \n\n")

    call <- match.call()
    txt <- scan(file,what="character",sep="\n",quiet=TRUE)
    txt <- gsub("\t"," ",txt)

                                        # read first infos
    info <- unlist(strsplit(txt[1],"([[:space:]]+)"))
                                        # npop <- as.numeric(info[1]) ## no longer used
    nloc <- as.numeric(info[2])

    loc.names <- txt[2:(nloc+1)]

                                        # build genotype matrix
    txt <- txt[-(1:(nloc+1))]
    txt <- .rmspaces(txt)
    txt <- sapply(1:length(txt),function(i) unlist(strsplit(txt[i],"([[:space:]]+)|([[:blank:]]+)")) )
    X <- t(txt)
    pop <- factor(X[,1])
    if(length(levels(pop)) == 1 ) pop <- NULL
    X <- X[,-1]

    colnames(X) <- loc.names
    rownames(X) <- 1:nrow(X)

    res <- df2genind(X=X,pop=pop,missing=missing, ploidy=2)
                                        # beware : fstat files do not yield ind names
    res@ind.names <- rep("",length(res@ind.names))
    names(res@ind.names) <- rownames(res@tab)
    res@call <- call

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.fstat





##########################
# Function read.genepop
##########################
read.genepop <- function(file,missing=NA,quiet=FALSE){
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(toupper(.readExt(file)) != "GEN") stop("File extension .gen expected")

    if(!quiet) cat("\n Converting data from a Genepop .gen file to a genind object... \n\n")

    prevcall <- match.call()

    txt <- scan(file,sep="\n",what="character",quiet=TRUE)
    if(!quiet) cat("\nFile description: ",txt[1], "\n")
    txt <- txt[-1]
    txt <- gsub("\t", " ", txt)

  # two cases for locus names:
  # 1) all on the same row, separated by ","
  # 2) one per row
  # ! spaces and tab allowed
  # a bug was reported by S. Devillard, occuring
  # when the two cases occur together,
  # that is:
  # loc1,
  # loc2,
  # ...

  ### former version
  #1
  #if(length(grep(",",txt[1])) > 0){
  #  loc.names <- unlist(strsplit(txt[1],","))
  #  loc.names <- gsub("^([[:blank:]]*)([[:space:]]*)","",loc.names)
  #  loc.names <- gsub("([[:blank:]]*)([[:space:]]*)$","",loc.names)
  #  nloc <- length(loc.names)

  #  txt <- txt[-1]
  #} else { #2
  #  nloc <- min(grep("POP",toupper(txt)))-1
  #  loc.names <- txt[1:nloc]
  #  loc.names <- gsub("^([[:blank:]]*)([[:space:]]*)","",loc.names)
  #  loc.names <- gsub("([[:blank:]]*)([[:space:]]*)$","",loc.names)

  #  txt <- txt[-(1:nloc)]
  #}

    ## new strategy (shorter): isolate the 'locus names' part and then parse it.
    locinfo.idx <- 1:(min(grep("POP",toupper(txt)))-1)
    locinfo <- txt[locinfo.idx]
    locinfo <- paste(locinfo,collapse=",")
    loc.names <- unlist(strsplit(locinfo,"([,]|[\n])+"))
    loc.names <- .rmspaces(loc.names)
    nloc <- length(loc.names)
    txt <- txt[-locinfo.idx]

    ## locus names have been retreived

    ## build the pop factor
    ## and correct the genotypes splited on more than 1 line
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))
    npop <- length(pop.idx)
    ## correction for splited genotype
    ## isolated by the absence of comma on a line not containing "pop"
    nocomma <- which(! (1:length(txt)) %in% grep(",",txt))
    splited <- nocomma[which(! nocomma %in% pop.idx)]
    if(length(splited)>0){
        for(i in sort(splited,dec=TRUE)){
            txt[i-1] <- paste(txt[i-1],txt[i],sep=" ")
        }
        txt <- txt[-splited]
    }
    ## end correction

    ## reevaluate pop index
    pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))

    txt[length(txt)+1] <- "POP"
    nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt)))-1
    pop <- factor(rep(1:npop,nind.bypop))

    txt <- txt[-c(pop.idx,length(txt))]

    temp <- sapply(1:length(txt),function(i) strsplit(txt[i],","))
    ## temp is a list with nind elements, first being ind. name and 2nd, genotype

    ind.names <- sapply(temp,function(e) e[1])
    ind.names <- .rmspaces(ind.names)
    ## individuals' name are now clean

    vec.genot <- sapply(temp,function(e) e[2])
    vec.genot <- .rmspaces(vec.genot)

    ## X is a individual x locus genotypes matrix
    X <- matrix(unlist(strsplit(vec.genot,"[[:space:]]+")),ncol=nloc,byrow=TRUE)

    rownames(X) <- ind.names
    colnames(X) <- loc.names

 ##  # correct X to fulfill the genetix format
##   f1 <- function(char){
##     paste("00", substr(char,1,1), "00", substr(char,2,2), sep="")
##   }

##   f2 <- function(char){
##     paste("0", substr(char,1,2), "0", substr(char,3,4), sep="")
##   }

##   if(all(nchar(X)==2)) {X <- apply(X,c(1,2),f1)}
##   if(all(nchar(X)==4)) {X <- apply(X,c(1,2),f2)}

    ## give right pop names
    ## beware: genepop takes the name of the last individual of a sample as this sample's name
    pop.names.idx <- cumsum(table(pop))
    pop.names <- ind.names[pop.names.idx]
    levels(pop) <- pop.names

    res <- df2genind(X=X,pop=pop,missing=missing, ploidy=2)
    res@call <- prevcall

    if(!quiet) cat("\n...done.\n\n")

    return(res)

} # end read.genepop





############################
# Function read.structure
############################
read.structure <- function(file, n.ind=NULL, n.loc=NULL,  onerowperind=NULL, col.lab=NULL, col.pop=NULL, col.others=NULL, row.marknames=NULL, NA.char="-9", pop=NULL, missing=NA, ask=TRUE, quiet=FALSE){

    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    if(!toupper(.readExt(file)) %in% c("STR","STRU")) stop("File extension .stru expected")

    ## set defaults for non optional arguments without default values
    if(!ask){
        if(is.null(col.lab)) col.lab <- as.integer(0)
        if(is.null(col.pop)) col.pop <- as.integer(0)
        if(is.null(row.marknames)) row.marknames <- as.integer(0)
    }

    ## required questions
    if(is.null(n.ind)){
        cat("\n How many genotypes are there? ")
        n.ind <- as.integer(readLines(n = 1))
    }

    if(is.null(n.loc)){
        cat("\n How many markers are there? ")
        n.loc <- as.integer(readLines(n = 1))
    }

    if(is.null(col.lab)){
        cat("\n Which column contains labels for genotypes ('0' if absent)? ")
        col.lab <- as.integer(readLines(n = 1))
    }

    if(is.null(col.pop)){
        cat("\n Which column contains the population factor ('0' if absent)? ")
        col.pop <- as.integer(readLines(n = 1))
    }

    if(is.null(col.others) & ask){
        cat("\n Which other optional columns should be read (press 'return' when done)? ")
        col.others <- scan(quiet=TRUE)
        if(length(col.others) == 0)  col.others <- NULL
    }

    if(is.null(row.marknames)){
        cat("\n Which row contains the marker names ('0' if absent)? ")
        row.marknames <- as.integer(readLines(n = 1))
    }

    if(is.null(onerowperind)){
        cat("\n Are genotypes coded by a single row (y/n)? ")
        onerowperind <- toupper(readLines(n = 1))
        if(onerowperind == "Y") {
            onerowperind <- TRUE
        } else {
            onerowperind <- FALSE
        }
    }

    if(is.null(NA.char)){
        cat("\n What is the code for missing data (default is '-9')? ")
        NA.char <- as.character(readLines(n = 1))
    }

    ## message to console
    if(!quiet) cat("\n Converting data from a STRUCTURE .stru file to a genind object... \n\n")

    ## read the file
    txt <- scan(file,sep="\n",what="character",quiet=TRUE)

    ## remove empty lines and spaces/tabs at the end of a line
    temp <- grep("^[[:space:]]*$",txt)
    if(length(temp) > 0) {
        txt <- txt[-temp]
    }

    txt <- gsub("([[:blank:]]+)$","",txt)

    ## isolate each useful component of the file
    ## matrix of data
    if(onerowperind) {
        n <- n.ind
        p <- 2*n.loc
    } else{
        n <- 2*n.ind
        p <- n.loc
    }

    lastline <- length(txt)
    mat <- txt[(lastline-n+1):lastline]
    mat <- t(as.data.frame(strsplit(mat,"[[:blank:]]+")))
    rownames(mat) <- 1:n
    gen <- mat[, (ncol(mat)-p+1):ncol(mat)]


    ## markers names
    if(row.marknames != 0) {
        loc.names <- .rmspaces(txt[row.marknames])
        loc.names <- unlist(strsplit(loc.names,"[[:blank:]]+"))
    } else {
        loc.names <- .genlab("L",n.loc)
    }

    ## genotypes labels
    if(col.lab !=0) {
        ind.names <- mat[, col.lab]
    } else {
        ind.names <- .genlab("",n.ind)
    }

    ## population factor
    if(col.pop !=0) {
        pop <- factor(mat[, col.pop])
    } else {
        pop <- NULL
    }

    ## other variables
    if(!is.null(col.others)){
        X.other <- mat[,col.others]
    }

    ## transformations if onerowperind is FALSE
    if(!onerowperind) {
        temp <- seq(1,n,by=2)
        ind.names <- ind.names[temp]
        if(length(ind.names) < n.ind) warning("Duplicated identifier for genotypes")
        pop <- pop[temp]
        if(exists("X.other")) X.other <- X.other[temp]

        ## make sur that all strings in gen have the same number of characters
        ncode <- max(nchar(gen))
        keepCheck <- any(nchar(gen) < ncode)

        while(keepCheck){
            mat0 <- matrix("", ncol=ncol(gen), nrow=nrow(gen))
            mat0[nchar(gen) < ncode] <- "0"
            gen <-  matrix(paste(mat0, gen, sep=""), nrow=nrow(mat0))
            keepCheck <- any(nchar(gen) < ncode)
        }

        ## reorder matrix of genotypes
        X <- t(sapply(temp, function(i) paste(gen[i,],gen[i+1,],sep="") ))

    } else { # else of "if(!onerowperind)"
        temp <- seq(1,p-1,by=2)
        X <- paste(gen[,temp] , gen[,temp+1], sep="")
        X <- matrix(X, nrow=n.ind)
    }

    ## replace missing values by NAs
    X <- gsub(NA.char,NA,X)
    rownames(X) <- ind.names
    colnames(X) <- loc.names

    res <- df2genind(X=X,pop=pop,missing=missing, ploidy=2)

    res@call <- match.call()

    if(exists("X.other")) {res@other <- list(X=X.other)}

    return(res)

}




#########################
# Function import2genind
#########################
import2genind <- function(file,missing=NA,quiet=FALSE, ...){
    ## if(!file.exists(file)) stop("Specified file does not exist.") <- not needed
    ext <- .readExt(file)
    ext <- toupper(ext)

    if(ext == "GTX")
        return(read.genetix(file,missing=missing,quiet=quiet))

    if(ext == "DAT")
        return(read.fstat(file,missing=missing,quiet=quiet))

    if(ext == "GEN")
        return(read.genepop(file,missing=missing,quiet=quiet))

    if(ext %in% c("STR","STRU"))
        return(read.structure(file,missing=missing,quiet=quiet, ...))

    ## evaluated only if extension is not supported
    cat("\n File format (",ext,") not supported.\n")
    cat("\nSupported formats are:\nGENETIX (.gtx) \nFSTAT (.dat) \nGenepop (.gen)\n \nSTRUCTURE (.str)\n")

    return(invisible())
}

