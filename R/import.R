##################################################################
# Fonctions designed to import files from other softwares
# into genind objects
#
# currently supported formats are :
# .gtx (GENETIX)
# .dat (Fstat)
# .gen (Genepop)
#
# Thibaut Jombart, avril 2006
# jombart@biomserv.univ-lyon1.fr
# 
##################################################################

#######################
# Function rmspaces
#######################
# removes spaces and tab at the begining and the end of each element of charvec
.rmspaces <- function(charvec){
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)","",charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$","",charvec)
    return(charvec)
}



########################################
# Function genetix2genind
# code based on previous ade4 functions
########################################
genetix2genind <- function(file=NULL,X=NULL,pop=NULL,missing=NA,quiet=FALSE) {
    if(is.null(X) == is.null(file)) stop("Either X or file must be provided")
    if(!quiet) cat("\n Converting data from GENETIX to a genind object... \n")

    if(!is.null(X)){
      if(is.data.frame(X)) X <- as.matrix(X)
      if (!inherits(X, "matrix")) stop ("X is not a matrix")
    }
    
    # eventually read from file
    if(!is.null(file)){
      if(!file.exists(file)) stop("Specified file does not exist.")
      if(toupper(substr(file,nchar(file)-2,nchar(file))) != "GTX") stop("File extension .gtx expected")
      # retrieve first infos
      nloc <- as.numeric(scan(file,nlines=1,what="character",quiet=TRUE)[1])
      npop <- as.numeric(scan(file,nlines=1,skip=1,what="character",quiet=TRUE)[1])
      txt <- scan(file,skip=2,what="character",sep="\n",quiet=TRUE)
      txt <- gsub("\t"," ",txt)
      loc.names <- txt[seq(1,by=2,length=nloc)]
      txt <- txt[-(1:(nloc*2))]

      # retrieve populations infos
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
      
      # retrieve genotypes infos
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
      
      # make a factor "pop" if there is more than one population
      pop <- factor(rep(pop.names,pop.nind))
    } # end if(!is.null(file))

    # now X exists whether file or X was provided by user

    # make the result
    res <- list()

    n <- nrow(X)
      
    # pop optionnelle
    if(!is.null(pop)){
      if(length(pop)!= n) stop("length of factor pop differs from nrow(X)")
      pop <- as.factor(pop)
    }
 
    # fonction pour corriger les cas a 4 caracteres, eliminer les autres, remplacer les NA
    checkcar <- function(cha) {
        n <- nchar(cha)
        if(n==6) return(cha)
        if(is.na(cha) || cha=="0") return("000000")
        if(n>6) stop(paste("data with more than 6 characters (",cha,")",sep=""))
        if(n==4) return(paste("0",substr(cha,1,2),"0",substr(cha,3,4),sep=""))
        if(n==2) return(paste("00",substr(cha,1,1),"00",substr(cha,2,2),sep=""))
        if(n<6) stop(paste("\n",cha,"is not interpretable"))  
     } # end checkcar
    
    if(any(nchar(X)) != 6) {X <- apply(X,c(1,2),checkcar)}
    # X contient a present seulement des données valides avec 6 caracteres

    # Erase entierely non-typed loci
    temp <- apply(X,2,function(c) all(c=="000000"))
    X <- X[,!temp]
    nloc <- ncol(X)
    loc.names <- colnames(X)
    
    # Erase entierely non-type individuals
    temp <- apply(X,1,function(c) all(c=="000000"))
    X <- X[!temp,]
    pop <- pop[!temp]
    n <- nrow(X)
    ind.names <- rownames(X)
    # note: if X is kept as a matrix, duplicate row names are no problem
    
    enumallel <- function (x) {
      w <- as.character(x)
      w1 <- substr(w,1,3)
      w2 <- substr(w,4,6)
      w3 <- sort(unique (c(w1,w2)))
      return(w3[w3!="000"])
    }

    loc.all <- lapply(1:ncol(X),function(i) enumallel(X[,i]))
    names(loc.all) <- loc.names
    # loc.all est une liste dont chaque element est un vecteur des alleles (ordonnes) d'un marqueur
    temp <- lapply(1:nloc, function(i) matrix(0,nrow=n,ncol=length(loc.all[[i]]),
       dimnames=list(NULL,loc.all[[i]])) )
    # note: keep rownames as NULL in case of duplicates
    
    names(temp) <- loc.names
    # temp est une liste dont chaque element est un tableau-marqueur (indiv x alleles)

    # remplissage des tableaux
    findall <- function(cha,loc.all){
      if(cha=="000") return(NULL)
      return(which(cha==loc.all))
    }

    for(i in 1:n){
      for(j in 1:nloc){
        all1pos <- findall(substr(X[i,j],1,3),loc.all[[j]])
        temp[[j]][i,all1pos] <- temp[[j]][i,all1pos] + 0.5
        all2pos <- findall(substr(X[i,j],4,6),loc.all[[j]])
        temp[[j]][i,all2pos] <- temp[[j]][i,all2pos] + 0.5
        if(is.null(c(all1pos,all2pos))) {temp[[j]][i,] <- NA}
      }
    }

    # beware: colnames are wrong when there is only one allele in a locus
    # right colnames are first generated
    nall <- unlist(lapply(temp,ncol))
    loc.rep <- rep(names(nall),nall)
    col.lab <- paste(loc.rep,unlist(loc.all,use.names=FALSE),sep=".")

    mat <- as.matrix(cbind.data.frame(temp))
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

    res <- as.genind( tab=mat, pop=pop, prevcall=prevcall )
    
    if(!quiet) cat("\n...done.\n\n")

    return(res)
} # end genetix2genind



##########################
# Function fstat2genind
##########################
fstat2genind <- function(file,missing=NA,quiet=FALSE){
  if(!file.exists(file)) stop("Specified file does not exist.")
  if(toupper(substr(file,nchar(file)-2,nchar(file))) != "DAT") stop("File extension .dat expected")

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

  res <- genetix2genind(X=X,pop=pop,missing=missing,quiet=TRUE)
  # beware : fstat files do not yield ind names
  res$ind.names <- rep("",length(res$ind.names))
  names(res$ind.names) <- rownames(res$tab)
  res$call <- call
  
  if(!quiet) cat("\n...done.\n\n")

  return(res)
  
} # end fstat2genind



##########################
# Function genepop2genind 
##########################
genepop2genind <- function(file,missing=NA,quiet=FALSE){
  if(!file.exists(file)) stop("Specified file does not exist.")
  if(toupper(substr(file,nchar(file)-2,nchar(file))) != "GEN") stop("File extension .gen expected")

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
  
  #1
  if(length(grep(",",txt[1])) > 0){
    loc.names <- unlist(strsplit(txt[1],","))
    loc.names <- gsub("^([[:blank:]]*)([[:space:]]*)","",loc.names)
    loc.names <- gsub("([[:blank:]]*)([[:space:]]*)$","",loc.names)
    nloc <- length(loc.names)

    txt <- txt[-1]
  } else { #2
    nloc <- min(grep("POP",toupper(txt)))-1
    loc.names <- txt[1:nloc]
    loc.names <- gsub("^([[:blank:]]*)([[:space:]]*)","",loc.names)
    loc.names <- gsub("([[:blank:]]*)([[:space:]]*)$","",loc.names)
  
    txt <- txt[-(1:nloc)]
  }

  # locus names have been retreived  

  # build the pop factor
  # and correct the genotypes splited on more than 1 line
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))
  npop <- length(pop.idx)
  # correction for splited genotype
  # isolated by the absence of comma on a line not containing "pop"
  nocomma <- which(! (1:length(txt)) %in% grep(",",txt))
  splited <- nocomma[which(! nocomma %in% pop.idx)]
  if(length(splited)>0){
    for(i in sort(splited,dec=TRUE)){
      txt[i-1] <- paste(txt[i-1],txt[i],sep=" ")
    }
    txt <- txt[-splited]
  }
  # end correction

  # reevaluate pop index
  pop.idx <- grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt))
  
  txt[length(txt)+1] <- "POP"
  nind.bypop <- diff(grep("^([[:space:]]*)POP([[:space:]]*)$",toupper(txt)))-1
  pop <- factor(rep(1:npop,nind.bypop))
  
  txt <- txt[-c(pop.idx,length(txt))]

  temp <- sapply(1:length(txt),function(i) strsplit(txt[i],","))
  # temp is a list with nind elements, first being ind. name and 2nd, genotype

  ind.names <- sapply(temp,function(e) e[1])
  ind.names <- .rmspaces(ind.names)
  # individuals' name are now clean

  vec.genot <- sapply(temp,function(e) e[2])
  vec.genot <- .rmspaces(vec.genot)
  
  # X is a individual x locus genotypes matrix
  X <- matrix(unlist(strsplit(vec.genot,"[[:space:]]+")),ncol=nloc,byrow=TRUE)
 
  rownames(X) <- ind.names
  colnames(X) <- loc.names

  # correct X to fulfill the genetix format
  f1 <- function(char){
    paste("00", substr(char,1,1), "00", substr(char,2,2), sep="")
  }

  f2 <- function(char){
    paste("0", substr(char,1,2), "0", substr(char,3,4), sep="")
  }

  if(all(nchar(X)==2)) {X <- apply(X,c(1,2),f1)}
  if(all(nchar(X)==4)) {X <- apply(X,c(1,2),f2)}

  # give right pop names
  # beware: genepop takes the name of the last individual of a sample as this sample's name
  pop.names.idx <- cumsum(table(pop))
  pop.names <- ind.names[pop.names.idx]
  levels(pop) <- pop.names
  
  res <- genetix2genind(X=X,pop=pop,missing=missing,quiet=TRUE)
  res$call <- prevcall
  
  if(!quiet) cat("\n...done.\n\n")

  return(res)
    
} # end genepop2genind



#########################
# Function import2genind
#########################
import2genind <- function(file,missing=NA,quiet=FALSE){
  if(!file.exists(file)) stop("Specified file does not exist.")
  ext <- toupper(substr(file,nchar(file)-2,nchar(file)))
  
  if(ext == "GTX")
    return(genetix2genind(file,missing=missing,quiet=quiet))

  if(ext == "DAT")
    return(fstat2genind(file,missing=missing,quiet=quiet))

  if(ext == "GEN")
    return(genepop2genind(file,missing=missing,quiet=quiet))

  # evaluated only if extension is not supported
  cat("\n File format (",ext,") not supported.\n")
  cat("\nSupported formats are:\nGENETIX (.gtx) \nFSTAT (.dat) \nGenepop (.gen)\n")
       
  return(invisible())    
}


