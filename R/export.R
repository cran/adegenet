############################################
#
# Functions to transform a genind object
# into other R classes
#
# Thibaut Jombart
# t.jombart@imperial.ac.uk
#
############################################



###########################
# Function genind2genotype
###########################
genind2genotype <- function(x,pop=NULL,res.type=c("matrix","list")){

  if(!is.genind(x)) stop("x is not a valid genind object")
  if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
  checkType(x)

  if(!require(genetics)) stop("genetics package is not required but not installed.")
  if(is.null(pop)) pop <- x@pop
  if(is.null(pop)) pop <- as.factor(rep("P1",nrow(x@tab)))
  res.type <- tolower(res.type[1])

  # make one table by locus from x@tab
  kX <- seploc(x,res.type="matrix")
  # kX is a list of nloc tables

  # function to recode a genotype in form "A1/A2" from frequencies
  recod <- function(vec,lab){
    if(all(is.na(vec))) return(NA)
    if(round(sum(vec),10) != 1) return(NA)
    temp <- c(which(vec==0.5),which(vec==1))
    if(length(temp)==0) return(NA)
    lab <- lab[temp]
    res <- paste(lab[1],lab[length(lab)],sep="/")
    return(res)
  }

  # function which converts data of a locus into a list of genotypes per population ## no longer used
  # f1 <- function(X){
  #  tapply(X,pop,function(mat) apply(mat,1,recod))
  #}

  # kGen is a list of nloc vectors of genotypes
  kGen <- lapply(1:length(kX), function(i) apply(kX[[i]],1,recod,x@all.names[[i]]))
  names(kGen) <- x@loc.names

  if(res.type=="list"){ # list type
    # each genotype is splited per population

    # correction of an error due to a change in as.genotype
    # error occurs in list type when a population is entierly untyped for a locus,
    # that is, all values are NA.
    res <- lapply(kGen,split,pop)

    f2 <- function(x){# x is a vector of character to be converted into genotype
      if(all(is.na(x))) return(NA)
      return(as.genotype(x))
    }
    res <- lapply(res,function(e) lapply(e,f2))
  } else if(res.type=="matrix"){ # matrix type
    res <- cbind.data.frame(kGen)
    res <- makeGenotypes(res,convert=1:ncol(res))
  } else stop("Unknown res.type requested.")

  return(res)
}





############################
# Function genind2hierfstat
############################
genind2hierfstat <- function(x,pop=NULL){
    ##  if(!inherits(x,"genind")) stop("x must be a genind object (see ?genind)")
    ##   invisible(validObject(x))
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(x@ploidy != as.integer(2)) stop("not implemented for non-diploid genotypes")
    checkType(x)

    if(is.null(pop)) pop <- x@pop
    if(is.null(pop)) pop <- as.factor(rep("P1",nrow(x@tab)))

    ## make one table by locus from x@tab
    kX <- seploc(x,res.type="matrix")
    ## kX is a list of nloc tables

    ## prepare allele names
    all.names <- x@all.names

    ## check the number of first 0 to remove from all.names
    nfirstzero <- attr(regexpr("^0*",unlist(all.names)),"match.length")
    nrmzero <- min(nfirstzero)

    for(i in 1:nrmzero) {
        all.names <- lapply(all.names,function(e) gsub("^0","",e))
    }

    ## function to recode a genotype in form "A1A2" (as integers) from frequencies
    recod <- function(vec,lab){
        if(all(is.na(vec))) return(NA)
        if(sum(vec) < 0) return(NA)
        temp <- which(vec!=0)
        lab <- lab[temp]
        res <- as.integer(paste(lab[1],lab[length(lab)],sep=""))
        return(res)
    }

                                        # kGen is a list of nloc vectors of genotypes
    kGen <- lapply(1:length(kX), function(i) apply(kX[[i]],1,recod,all.names[[i]]))
    res <- cbind(as.numeric(pop),as.data.frame(kGen))
    colnames(res) <- c("pop",x@loc.names)

    return(res)
}




#####################
# Function genind2df
#####################
genind2df <- function(x, pop=NULL, sep="", usepop=TRUE, oneColPerAll=FALSE){

  if(!is.genind(x)) stop("x is not a valid genind object")
  ## checkType(x)

  if(is.null(pop)) {
      pop <- x@pop
      levels(pop) <- x@pop.names
  }

  if(oneColPerAll){
      sep <- "/"
  }

  ## PA case ##
  if(x@type=="PA"){
      temp <- truenames(x)
      if(is.list(temp) & usepop){
          res <- cbind.data.frame(pop=temp[[2]],temp[[1]])
      } else{
          if(is.list(temp)) {
              res <- temp[[1]]
          } else{
              res <- temp
          }
      }

      return(res) # exit here
  }

  ## codom case ##
  # make one table by locus from x@tab
  kX <- seploc(x,res.type="matrix")
  kX <- lapply(kX, function(X) round(X*x@ploidy)) # take data as numbers of alleles
  ## (kX is a list of nloc tables)

  ## function to recode a genotype in form "A1[sep]...[sep]Ak" from frequencies
  recod <- function(vec,lab){
      if(any(is.na(vec))) return(NA)
      res <- paste( rep(lab,vec), collapse=sep)
      return(res)
  }


  # kGen is a list of nloc vectors of genotypes
  kGen <- lapply(1:length(kX), function(i) apply(kX[[i]],1,recod,x@all.names[[i]]))
  names(kGen) <- x@loc.names

  ## if use one column per allele
  if(oneColPerAll){
      f1 <- function(vec){ # to repeat NA with seperators
          vec[is.na(vec)] <- paste(rep("NA",x@ploidy), collapse=sep)
          return(vec)
      }
      temp <- lapply(kGen, f1)
      temp <- lapply(temp, strsplit,sep)

      res <- lapply(temp, function(e) matrix(unlist(e), ncol=x@ploidy, byrow=TRUE))
      res <- data.frame(res,stringsAsFactors=FALSE)
      names(res) <- paste(rep(locNames(x),each=x@ploidy), 1:x@ploidy, sep=".")

      ## handle pop here
      if(!is.null(pop) & usepop) res <- cbind.data.frame(pop,res,stringsAsFactors=FALSE)

      return(res) # exit here
  } # end if oneColPerAll

  ## build the final data.frame
  res <- cbind.data.frame(kGen,stringsAsFactors=FALSE)

  ## handle pop here
  if(!is.null(pop) & usepop) res <- cbind.data.frame(pop,res,stringsAsFactors=FALSE)

  return(res)
}
