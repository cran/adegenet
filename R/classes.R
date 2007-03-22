########################################################################
# Ensemble de fonctions pour le traitement des données génétiques
# dans le cadre de l'interface avec des méthodes d'analyse multivariées
# du package ade4 de R.
#
# Thibaut Jombart, Février 2007
# jombart@biomserv.univ-lyon1.fr
########################################################################

###############################
# Two classes of R object are 
# defined :
# genind - genotypes of individuals
# genpop - allelic frequencies of populations
###############################


#############################################################################################
#############################################################################################
# AUXILIARY FUNCTIONS
#############################################################################################
#############################################################################################



#######################
# Function rmspaces
#######################
# removes spaces and tab at the begining and the end of each element of charvec
.rmspaces <- function(charvec){
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)","",charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$","",charvec)
    return(charvec)
}



#######################################
# Function .is.gen
# used in both is.genind and is.genpop
#######################################
.is.gen <- function(x){
  if(!is.list(x)) {
    cat("\nnot a list\n")
    return(FALSE)
  }
  
# tests sur l'existence des éléments
  if( is.null(x$tab) ||
     is.null(x$loc.names) ||
     is.null(x$loc.fac) ||
     is.null(x$loc.nall) ||
     is.null(x$all.names) ||
     is.null(x$call)) {
       cat("\nMissing elements.\n")
       return(FALSE)
     }

# tests sur la nature des éléments
  if(is.data.frame(x$tab)) x$tab <- as.matrix(x$tab)
  if(!is.matrix(x$tab)) {
    cat("\ntab not a matrix\n")
    return(FALSE)
  }
  
  if(!is.character(x$loc.names)) {
    cat("\nloc.names not a character vector\n")
    return(FALSE)
  }
  
  if(!is.factor(x$loc.fac)) {
    cat("\nloc.fac not a factor\n")
    return(FALSE)
  }
  
  if(!is.numeric(x$loc.nall)) {
    cat("\nloc.nall not a numeric vector\n")
    return(FALSE)
  }
  
  if(!is.list(x$all.names)) {
    cat("\nall.names is not a list\n")
    return(FALSE)
  }

  if(!all(lapply(x$all.names,is.character))) {
    cat("\nall.names does not contains only character vectors\n")
    return(FALSE)
  }

  # tests des dimensions des éléments 
  n <- nrow(x$tab)
  p <- ncol(x$tab)
  k <- length(unique(x$loc.names))

  if(length(x$loc.fac) != p) {
    cat("\ninvalid length for loc.fac\n")
    return(FALSE)
  }

  if(length(levels(x$loc.fac)) != k) {
    cat("\ninvalid number of levels in loc.fac\n")
    return(FALSE)
  }

  if(length(x$loc.nall) != k) {
    cat("\ninvalid length in loc.nall\n")
    return(FALSE)
  }

  if(length(unlist(x$all.names)) != p) {
    cat("\ninvalid length in all.names\n")
    return(FALSE)
  }

  return(TRUE)

}# end .is.gen



###################
# Function .genlab
###################
# recursive function to have labels of constant length
# base = a character string
# n = number of labels
.genlab <- function(base, n) {
  f1 <- function(cha,n){
    if(nchar(cha)<n){
      cha <- paste("0",cha,sep="")
      return(f1(cha,n))
    } else {return(cha)}
  }
  w <- as.character(1:n)
  max0 <- max(nchar(w))
  w <- sapply(w, function(cha) f1(cha,max0))
  return(paste(base,w,sep=""))
}




#############################################################################################
#############################################################################################
# CLASS FUNCTIONS
#############################################################################################
#############################################################################################



#####################
# Function .is.genind
#####################
is.genind <- function(x){
  if(!inherits(x,"genind")) return(FALSE)

  if(!.is.gen(x)) return(FALSE)

  if(!is.null(x$pop)){
    if(is.character(x$pop) || is.numeric(x$pop)) x$pop <- as.factor(x$pop)

    if(!is.factor(x$pop)) {
      cat("\npop is given but not a factor\n")
      return(FALSE)
    }

    if(length(x$pop) != nrow(x$tab)) {
      cat("\npop is given but invalid length\n")
      return(FALSE)
    }
  }

  if(is.null(x$ind.names)) {
    cat("\nind.names is missing\n")
    return(FALSE)
  }

  if(!is.character(x$ind.names)) {
      cat("\nind.names is not a character vector\n")
      return(FALSE)
    }

   if(length(x$ind.names) != nrow(x$tab)) {
      cat("\ninvalid length in ind.names\n")
      return(FALSE)
    }

  if(!is.null(x$pop.names)){
    if(is.null(x$pop)) warning("\npop.names is given but pop is not\n")
    if(!is.character(x$pop.names)) {
      cat("\npop.names is given but not a character vector\n")
      return(FALSE)
    }
    if(length(x$pop.names) != length(levels(x$pop))) {
      cat("\npop.names is given but invalid length\n")
      return(FALSE)
    }
  }

  return(TRUE)
  
}# end is.genind



#####################
# Function is.genpop
#####################
is.genpop <- function(x){
  if(!inherits(x,"genpop")) return(FALSE)

  if(!.is.gen(x)) return(FALSE)
  
  if(is.null(x$pop.names)) {
    cat("\npop.names missing\n")
    return(FALSE)
  }

  if(!is.character(x$pop.names)) {
    cat("\npop.names is not a character vector\n")
    return(FALSE)
  }

  if(length(x$pop.names) != nrow(x$tab)) {
    cat("\ninvalid length in pop.names\n")
    return(FALSE)
  }
  
  return(TRUE)
  
}# end is.genpop



#####################
# Function as.genind
#####################
as.genind <- function(tab=NULL,pop=NULL,prevcall=NULL){

  if(is.null(tab)) stop("tab not provided.")
  X <- as.matrix(tab)
  if(is.null(colnames(X))) stop("tab columns have no name.")
  if(is.null(rownames(X))) {rownames(X) <- 1:nrow(X)}
 
  
  # labels for individuals
  nind <- nrow(X)
  ind.names <- .rmspaces(rownames(X))
  ind.codes <- .genlab("", nind)
  names(ind.names) <- ind.codes
  
  # labels for loci
  # and loc.nall
  temp <- colnames(X)
  temp <- gsub("[.].*$","",temp)
  temp <- .rmspaces(temp)
  # beware !!! Function 'table' gives ordred output.
  loc.names <- unique(temp)
  loc.nall <-  table(temp)[match(loc.names,names(table(temp)))]    

  nloc <- length(loc.names)
  loc.codes <- .genlab("L",nloc)

  names(loc.names) <- loc.codes

  names(loc.nall) <- loc.codes

  # loc.fac
  loc.fac <- rep(loc.codes,loc.nall)

  # alleles name
  temp <- colnames(X)
  temp <- gsub("^.*[.]","",temp)
  temp <- .rmspaces(temp)
  all.names <- split(temp,loc.fac)
  all.codes <- lapply(all.names,function(e) .genlab("",length(e)))
  for(i in 1:length(all.names)){
    names(all.names[[i]]) <- all.codes[[i]]
  }
  
  rownames(X) <- ind.codes
  colnames(X) <- paste(loc.fac,unlist(all.codes),sep=".")
  loc.fac <- as.factor(loc.fac)
  
  res <- list( tab=X, ind.names=ind.names, loc.names=loc.names,
              loc.nall=loc.nall, loc.fac=loc.fac, all.names=all.names )
  
  # populations name (optional)
  if(!is.null(pop)) {
    pop.lab <- .genlab("P",length(levels(pop)) )
    temp <- pop
    levels(pop) <- pop.lab
    res$pop <- pop
    res$pop.names <- as.character(levels(temp))
    names(res$pop.names) <- as.character(levels(res$pop))
  }

  if(is.null(prevcall)) {prevcall <- match.call()}
  res$call <- prevcall

  class(res) <- "genind"
  
  return(res)
  
} # end as.genind



#####################
# Function as.genpop
#####################
as.genpop <- function(tab=NULL,prevcall=NULL){

  if(is.null(tab)) stop("tab not provided.")
  X <- as.matrix(tab)
  if(is.null(colnames(X))) stop("tab columns have no name.")
  if(is.null(rownames(X))) {rownames(X) <- 1:nrow(X)}
 
  # labels for populations
  npop <- nrow(X)
  pop.names <- .rmspaces(rownames(X))
  pop.codes <- .genlab("P", npop)
  names(pop.names) <- pop.codes
  
  # labels for loci
  # and loc.nall
  temp <- colnames(X)
  temp <- gsub("[.].*$","",temp)
  temp <- .rmspaces(temp)
  # beware !!! Function 'table' gives ordred output.
  loc.names <- unique(temp)
  loc.nall <-  table(temp)[match(loc.names,names(table(temp)))]    

  nloc <- length(loc.names)
  loc.codes <- .genlab("L",nloc)

  names(loc.names) <- loc.codes

  names(loc.nall) <- loc.codes

  # loc.fac
  loc.fac <- rep(loc.codes,loc.nall)

  # alleles name
  temp <- colnames(X)
  temp <- gsub("^.*[.]","",temp)
  temp <- .rmspaces(temp)
  all.names <- split(temp,loc.fac)
  all.codes <- lapply(all.names,function(e) .genlab("",length(e)))
  for(i in 1:length(all.names)){
    names(all.names[[i]]) <- all.codes[[i]]
  }
  
  rownames(X) <- pop.codes
  colnames(X) <- paste(loc.fac,unlist(all.codes),sep=".")
  loc.fac <- as.factor(loc.fac)
  
  res <- list( tab=X, pop.names=pop.names, loc.names=loc.names,
              loc.nall=loc.nall, loc.fac=loc.fac, all.names=all.names )
  
  if(is.null(prevcall)) {prevcall <- match.call()}
  res$call <- prevcall

  class(res) <- "genpop"
  
  return(res)
  
} # end as.genpop



########################
# Function print.genind
########################
print.genind <- function(x,...){
  cat("   #####################\n")
  cat("   ### Genind object ### \n")
  cat("   #####################")
  cat("\n- genotypes of individuals - \n")
  cat("\nclass: ")
  print(class(x))
  cat("\n$call: ")
  print(x$call)

  n <- nrow(x$tab)
  p <- ncol(x$tab)
  len <- 7

  cat("\n$tab: ", nrow(x$tab), "x", ncol(x$tab), "matrix of genotypes\n" )
  head(x$tab[,1:min(p,5)],3)

  cat("\n$ind.names: vector of ", length(x$ind.names), "individual names")
  cat("\n$loc.names: vector of ", length(x$loc.names), "locus names")
  cat("\n$loc.nall: number of alleles per locus")
  cat("\n$loc.fac: locus factor for the ", ncol(x$tab), "columns of $tab")
  cat("\n$all.names: list of ", length(x$all.names), "components yielding allele names for each locus")
  if(!is.null(x$pop)) {
    cat("\n\nOptionnal contents: ")
    cat("\n$pop: factor giving the population of each individual")
    len=len+1
  }
  if(!is.null(x$pop.names)) {
    cat("\n$pop.names: vector giving the names of the populations")
    len=len+1
  }
 
  cat("\n\nother elements: ")
  if (length(names(x)) > len) cat(paste("$",names(x)[(len+1):(length(names(x)))],sep=""), "\n")
  else cat("NULL\n")
  return(len)
} #end print.genind 



########################
# Function print.genpop
########################
print.genpop <- function(x,...){
  cat("       #####################\n")
  cat("       ### Genpop object ### \n")
  cat("       #####################")
  cat("\n- Alleles counts for populations - \n")
  cat("\nclass: ")
  print(class(x))
  cat("\n$call: ")
  print(x$call)

  n <- nrow(x$tab)
  p <- ncol(x$tab)

  cat("\n$tab: ", nrow(x$tab), "x", ncol(x$tab), "matrix of alleles counts\n" )
  head(x$tab[,1:min(p,5)],3)

  cat("\n$pop.names: vector of ", length(x$pop.names), "population names")
  cat("\n$loc.names: vector of ", length(x$loc.names), "locus names")
  cat("\n$loc.nall: number of alleles per locus")
  cat("\n$loc.fac: locus factor for the ", ncol(x$tab), "columns of $tab")
  cat("\n$all.names: list of ", length(x$all.names), "components yielding allele names for each locus")
 
  cat("\n\nother elements: ")
  if (length(names(x)) > 7) cat(paste("$",names(x)[8:(length(names(x)))],sep=""), "\n")
  else cat("NULL\n")
} # end print.genpop



##########################
# Function summary.genind
##########################
summary.genind <- function(object, ...){
  x <- object
  if(!is.genind(x)) stop("To be used with a genind object")
  if(is.null(x$pop)){
    x$pop <- rep(1,nrow(x$tab))
    x$pop.names <- ""
    names(x$pop.names) <- "P1"
  }

  res <- list()

  res$N <- nrow(x$tab)

  res$pop.eff <- as.numeric(table(x$pop))
  names(res$pop.eff) <- names(x$pop.names)

  res$loc.nall <- x$loc.nall

  temp <- genind2genpop(x,quiet=TRUE)$tab

  res$pop.nall <- apply(temp,1,function(r) sum(r!=0,na.rm=TRUE))

  res$NA.perc <- 100*sum(is.na(x$tab))/prod(dim(x$tab))

  # auxiliary function to compute observed heterozygosity 
  temp <- seploc(x,truenames=FALSE) 
  f1 <- function(tab){
    H <- sum(tab==0.5,na.rm=TRUE)/(2*nrow(tab))
    return(H)
  }

  res$Hobs <- unlist(lapply(temp,f1))

  # auxiliary function to compute expected heterozygosity
  # freq is a vector of frequencies
  f2 <- function(freq){
    H <- 1-sum(freq*freq,na.rm=TRUE)
    return(H)
  }

  temp <- genind2genpop(x,pop=rep(1,nrow(x$tab)),quiet=TRUE)
  temp <- makefreq(temp,quiet=TRUE)$tab
  temp.names <- colnames(temp)
  temp <- as.vector(temp)
  names(temp) <- temp.names
  temp <- split(temp,x$loc.fac)
  # temp is a list of alleles frequencies (one element per locus)

  res$Hexp <- unlist(lapply(temp,f2))  

  # print to screen
  listlab <- c("# Total number of genotypes: ",
               "# Population sample sizes: ",
               "# Number of alleles per locus: ",
               "# Number of alleles per population: ",
               "# Percentage of missing data: ",
               "# Observed heterozygosity: ",
               "# Expected heterozygosity: ")
  cat("\n",listlab[1],res[[1]],"\n")
  for(i in 2:7){
    cat("\n",listlab[i],"\n")
    print(res[[i]])
  }
  
  return(invisible(res))
} # end summary.genind



##########################
# Function summary.genpop
##########################
summary.genpop <- function(object,...){
  x <- object
  if(!is.genpop(x)) stop("To be used with a genpop object")
  res <- list()

  res$npop <- nrow(x$tab)

  res$loc.nall <- x$loc.nall

  res$pop.nall <- apply(x$tab,1,function(r) sum(r>0,na.rm=TRUE))

  res$NA.perc <- 100*sum(is.na(x$tab))/prod(dim(x$tab))

  # print to screen
  listlab <- c("# Number of populations: ",
               "# Number of alleles per locus: ",
               "# Number of alleles per population: ",
               "# Percentage of missing data: ")
  cat("\n",listlab[1],res[[1]],"\n")
  for(i in 2:4){
    cat("\n",listlab[i],"\n")
    print(res[[i]])
  }
  
  return(invisible(res))

} # end summary.genpop



#########################
# Function genind2genpop
#########################
genind2genpop <- function(x,pop=NULL,missing=NA,quiet=FALSE){
  if(!is.genind(x)) stop("x must be a genind object")
  if(is.null(x$pop) && is.null(pop)) stop("pop is not provided either in x or in pop")

  if(!quiet) cat("\n Converting data from a genind to a genpop object... \n")

  res <- list()
  
  # choose pop argument over x$pop
   if(!is.null(pop)) {
    if(length(pop) != nrow(x$tab)) stop("inconsistent length for factor pop")
    pop <- as.factor(pop)
  } else {
    pop <- x$pop
    if(!is.null(x$pop.names)) levels(pop) <- x$pop.names # restore real names
  }

  # make generic pop labels, store real pop names
  pop.names <- levels(pop)

  # tabcount is a matrix pop x alleles, counting alleles per pop
  # *2 to have alleles count
  f1 <- function(v){
    if(all(is.na(v))) return(NA) else return(sum(v,na.rm=TRUE))
  }

  f2 <- function(v){
    if(all(is.na(v)) || sum(v,na.rm=TRUE)==0) return(NA)
    return(v/(sum(v,na.rm=TRUE)))       
  }
  
  tabcount <- 2* apply(x$tab,2,function(c) tapply(c,pop,f1))
  # restitue matrix class when only one pop
  if(is.null(dim(tabcount))) {
    lab.col <- names(tabcount)
    tabcount <- matrix(tabcount,nrow=1)
    colnames(tabcount) <- lab.col
  }
  meancol <- apply(tabcount,2,function(c) mean(c,na.rm=TRUE))

  # NA treatment
  # Treatment when missing='REPLACE':
  # if allele 'j' of locus 'k' in pop 'i' is missing, replace the NA by a number 'x' so that
  # the frequency 'x/s' ('s' being the number of observations in 'k' ) equals the frequency 'f'
  # computed on the whole data (i.e. considering all pop as one)
  # Then x must verify:
  # x/s = f(1-f) => x=f(1-f)s
  #
  # - eff.pop is a pop x locus matrix giving the corresponding sum of observations (i.e., 's')
  # - temp is the same table but duplicated for all alleles
  # - odd.vec is the vector of 'f(1-f)'
  # - count.replace is a pop x alleles table yielding appropriate replacement numbers (i.e., 'x')

  if(!is.na(missing) && any(is.na(tabcount))){
    if(missing==0) tabcount[is.na(tabcount)] <- 0
    if(toupper(missing)=="REPLACE") {
    eff.pop <- t(apply(tabcount,1,function(r) tapply(r,x$loc.fac,sum,na.rm=TRUE)))
    temp <- t(apply(eff.pop,1,function(r) rep(r,table(x$loc.fac))))

    freq.allpop <- apply(tabcount,2,sum,na.rm=TRUE)
    freq.allpop <- unlist(tapply(freq.allpop,x$loc.fac,f2))
    odd.vec <- freq.allpop/(1-freq.allpop)
  
    count.replace <- t(apply(temp,1,function(r) r*odd.vec))

    tabcount[is.na(tabcount)] <- count.replace[is.na(tabcount)]
    }
  } # end of NA treatment

  # make final object
  temp <- paste(rep(x$loc.names,x$loc.nall),unlist(x$all.names),sep=".")
  colnames(tabcount) <- temp

  prevcall <- match.call()
  
  res <- as.genpop(tab=tabcount, prevcall=prevcall)

  if(!quiet) cat("\n...done.\n\n")

  return(res)
  
} # end genind2genpop

