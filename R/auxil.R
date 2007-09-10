###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################



#####################
# Function truenames
#####################
truenames <- function(x){

  if(is.genind(x)){
    return(.truenames.genind(x))
  }

  if(is.genpop(x)){
    return(.truenames.genpop(x))
  }

  stop("To be used with genind or genpop objects.")
}



#############################
# Function .truenames.genind
#############################
.truenames.genind <- function(x){

  X <- x$tab
  if(!all(x$ind.names=="")) {rownames(X) <- x$ind.names}

  labcol <- rep(x$loc.names,x$loc.nall)
  labcol <- paste(labcol,unlist(x$all.names),sep=".")
  colnames(X) <- labcol

  if(!is.null(x$pop)){
    pop <- x$pop
    levels(pop) <- x$pop.names
    return(list(tab=X,pop=pop))
  }

  return(X)
}



#############################
# Function .truenames.genpop
#############################
.truenames.genpop <- function(x){

  X <- x$tab
  if(!all(x$pop.names=="")) {rownames(X) <- x$pop.names}

  labcol <- rep(x$loc.names,x$loc.nall)
  labcol <- paste(labcol,unlist(x$all.names),sep=".")
  colnames(X) <- labcol

  return(X)
}



#####################
# Function seploc
#####################
seploc <- function(x,truenames=FALSE){

  if(is.genind(x)){
    return(.seploc.genind(x,truenames))
  }

  if(is.genpop(x)){
    return(.seploc.genpop(x,truenames))
  }

  stop("To be used with genind or genpop objects.")
}



##########################
# Function .seploc.genind
##########################
.seploc.genind <- function(x,truenames=FALSE){
  if(!is.genind(x)) stop("Wrong genind object.")
  
  temp <- x$loc.fac
  nloc <- length(levels(temp))
  levels(temp) <- 1:nloc

  kX <- list()
  
  for(i in 1:nloc){
    kX[[i]] <- matrix(x$tab[,temp==i],ncol=x$loc.nall[i])

    if(!truenames){
      rownames(kX[[i]]) <- rownames(x$tab)
      colnames(kX[[i]]) <- paste(names(x$loc.names)[i],names(x$all.names[[i]]),sep=".")
    }else{
      rownames(kX[[i]]) <- x$ind.names
      colnames(kX[[i]]) <- paste(x$loc.names[i],x$all.names[[i]],sep=".")
    }
  }

  if(truenames) {
    names(kX) <- x$loc.names
  } else{
    names(kX) <- names(x$loc.names)
  }

  return(kX)  
}



##########################
# Function .seploc.genpop
##########################
.seploc.genpop <- function(x,truenames=FALSE){
  if(!is.genpop(x)) stop("Wrong genpop object.")
  
  temp <- x$loc.fac
  nloc <- length(levels(temp))
  levels(temp) <- 1:nloc

  kX <- list()
  
  for(i in 1:nloc){
    kX[[i]] <- matrix(x$tab[,temp==i],ncol=x$loc.nall[i])

    if(!truenames){
      rownames(kX[[i]]) <- rownames(x$tab)
      colnames(kX[[i]]) <- paste(names(x$loc.names)[i],names(x$all.names[[i]]),sep=".")
    }else{
      rownames(kX[[i]]) <- x$pop.names
      colnames(kX[[i]]) <- paste(x$loc.names[i],x$all.names[[i]],sep=".")
    }
  }

  if(truenames) {
    names(kX) <- x$loc.names
  } else{
    names(kX) <- names(x$loc.names)
  }

  return(kX)  
}



#######################
# Function adegenetWeb
#######################
adegenetWeb <- function(){
  cat("Opening url \"http://pbil.univ-lyon1.fr/software/adegenet/\" ...\n")
  browseURL("http://pbil.univ-lyon1.fr/software/adegenet/")
}
