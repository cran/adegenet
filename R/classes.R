########################################################################
# adegenet classes definitions. All classes are S4.
#
# Thibaut Jombart, November 2007
# jombart@biomserv.univ-lyon1.fr
########################################################################

###############################
# Two classes of R object are 
# defined :
# gen - common part to genind and genpop
# genind - genotypes of individuals
# genpop - allelic frequencies of populations
###############################


###############################################################
###############################################################
# AUXILIARY FUNCTIONS
###############################################################
###############################################################



#######################
# Function rmspaces
#######################
# removes spaces and tab at the begining and the end of each element of charvec
.rmspaces <- function(charvec){
    charvec <- gsub("^([[:blank:]]*)([[:space:]]*)","",charvec)
    charvec <- gsub("([[:blank:]]*)([[:space:]]*)$","",charvec)
    return(charvec)
}



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



###############################################################
###############################################################
# CLASSES DEFINITION
###############################################################
###############################################################

#.initAdegenetClasses <- function(){


####################
# Unions of classes
####################
setClassUnion("listOrNULL", c("list","NULL"))
setClassUnion("factorOrNULL", c("factor","NULL"))
setClassUnion("charOrNULL", c("character","NULL"))
setClassUnion("callOrNULL", c("call","NULL"))
setClassUnion("intOrNum", c("integer","numeric"))



####################
# virtual class gen
####################
.gen.valid <- function(object){
  # this function tests only the consistency
  # of the length of each component
  p <- ncol(object@tab)
  k <- length(unique(object@loc.names))

  if(length(object@loc.fac) != p) {
    cat("\ninvalid length for loc.fac\n")
    return(FALSE)
  }

  if(length(levels(object@loc.fac)) != k) {
    cat("\ninvalid number of levels in loc.fac\n")
    return(FALSE)
  }

  if(length(object@loc.nall) != k) {
    cat("\ninvalid length in loc.nall\n")
    return(FALSE)
  }

  if(length(unlist(object@all.names)) != p) {
    cat("\ninvalid length in all.names\n")
    return(FALSE)
  }

  return(TRUE)

}# end .gen.valid


setClass("gen", representation(tab = "matrix",
                               loc.names = "character",
                               loc.fac = "factor",
                               loc.nall = "intOrNum",
                               all.names = "list",
                               call = "callOrNULL",
                               "VIRTUAL"),
         prototype(tab=matrix(ncol=0,nrow=0), loc.nall=integer(0), call=NULL))

setValidity("gen", .gen.valid)





########################
# virtual class indInfo
########################
setClass("indInfo", representation(ind.names = "character",
                                  pop = "factorOrNULL",
                                  pop.names = "charOrNULL",
                                  other = "listOrNULL", "VIRTUAL"),
         prototype(pop=NULL, pop.names = NULL, other = NULL))





###############
# Class genind
###############
.genind.valid <- function(object){
    if(!.gen.valid(object)) return(FALSE)
    
    if(length(object@ind.names) != nrow(object@tab)) {
        cat("\ninvalid length in ind.names\n")
        return(FALSE)
    }
    
    if(!is.null(object@pop)){ # check pop
        
        if(length(object@pop) != nrow(object@tab)) {
            cat("\npop is given but has invalid length\n")
            return(FALSE)
        }
        
        if(is.null(object@pop.names)) {
            cat("\npop is provided without pop.names")
        }  
        
        if(length(object@pop.names) != length(levels(object@pop))) {
            cat("\npop.names has invalid length\n")
            return(FALSE)
        }
    } # end check pop
    
    return(TRUE)
} #end .genind.valid

setClass("genind", contains=c("gen", "indInfo"))
setValidity("genind", .genind.valid)



########################
# virtual class popInfo
########################
setClass("popInfo", representation(pop.names = "character", other = "listOrNULL", "VIRTUAL"),
         prototype(other = NULL))



###############
# Class genpop
###############
.genpop.valid <- function(object){
    if(!.gen.valid(object)) return(FALSE)
    if(length(object@pop.names) != nrow(object@tab)) {
        cat("\ninvalid length in pop.names\n")
        return(FALSE)
    }
    
    return(TRUE)
} #end .genpop.valid

setClass("genpop", contains=c("gen", "popInfo"))
setValidity("genpop", .genpop.valid)







###############################################################
###############################################################
# MAIN CLASS METHODS
###############################################################
###############################################################



#################
# Function names
#################
setMethod("names", signature(x = "genind"), function(x){
  temp <- rev(names(attributes(x)))[-1]
  return(rev(temp))
})

setMethod("names", signature(x = "genpop"), function(x){
  temp <- rev(names(attributes(x)))[-1]
  return(rev(temp))
})





##################
# Function genind
##################
## constructor of a genind object
genind <- function(tab,pop=NULL,prevcall=NULL){

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
  loc.nall <- as.integer(loc.nall)

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
  
  # This was used in S3 version
  #
  #res <- list( tab=X, ind.names=ind.names, loc.names=loc.names,
  #            loc.nall=loc.nall, loc.fac=loc.fac, all.names=all.names )

  # Ideally I should use an 'initialize' method here
  res <- new("genind")
  res@tab <- X
  res@ind.names <- ind.names
  res@loc.names <- loc.names
  res@loc.nall <- loc.nall
  res@loc.fac <- loc.fac
  res@all.names <- all.names
  
  # populations name (optional)
  # beware, keep levels of pop sorted in
  # there order of appearance
  if(!is.null(pop)) {
      # convert pop to a factor if it is not
      if(!is.factor(pop)) {pop <- factor(pop)}
      pop.lab <- .genlab("P",length(levels(pop)) )
      # put pop levels in appearance order
      pop <- as.character(pop)
      pop <- factor(pop, levels=unique(pop))
      temp <- pop
      # now levels are correctly ordered
      levels(pop) <- pop.lab
      res@pop <- pop
      pop.names <- as.character(levels(temp))
      names(pop.names) <- as.character(levels(res@pop))
      res@pop.names <- pop.names
  }

  if(is.null(prevcall)) {prevcall <- match.call()}
  res@call <- prevcall
 
  return(res)
  
} # end genind

######################
# alias for as.genind
######################
as.genind <- genind



##################
# Function genpop
##################
genpop <- function(tab,prevcall=NULL){

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
  loc.nall <- as.integer(loc.nall)

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

  # Old S3 version
  #
  #res <- list( tab=X, pop.names=pop.names, loc.names=loc.names,
  #            loc.nall=loc.nall, loc.fac=loc.fac, all.names=all.names )

  res <- new("genpop")

  res@tab <- X
  res@pop.names <- pop.names
  res@loc.names <- loc.names
  res@loc.nall <- loc.nall
  res@loc.fac <- loc.fac
  res@all.names <- all.names
 
  if(is.null(prevcall)) {prevcall <- match.call()}
  res@call <- prevcall
  
  return(res)
  
} # end genpop



######################
# alias for as.genpop
######################
as.genpop <- genpop




##########################
# Method show for genind
##########################
setMethod ("show", "genind", function(object){
  x <- object
  cat("\n")
  cat("   #####################\n")
  cat("   ### Genind object ### \n")
  cat("   #####################")
  cat("\n- genotypes of individuals - \n")
  cat("\nS4 class: ", as.character(class(x)))
  
  cat("\n@call: ")
  print(x@call)

  p <- ncol(x@tab)
  len <- 7

  cat("\n@tab: ", nrow(x@tab), "x", ncol(x@tab), "matrix of genotypes\n" )
 
  cat("\n@ind.names: vector of ", length(x@ind.names), "individual names")
  cat("\n@loc.names: vector of ", length(x@loc.names), "locus names")
  cat("\n@loc.nall: number of alleles per locus")
  cat("\n@loc.fac: locus factor for the ", ncol(x@tab), "columns of @tab")
  cat("\n@all.names: list of ", length(x@all.names), "components yielding allele names for each locus")

  cat("\n\nOptionnal contents: ")
  cat("\n@pop: ", ifelse(is.null(x@pop), "- empty -", "factor giving the population of each individual"))
  cat("\n@pop.names: ", ifelse(is.null(x@pop.names), "- empty -", "factor giving the population of each individual"))
 
  cat("\n\n@other: ")
  if(!is.null(x@other)){
    cat("a list containing: ")
    cat(ifelse(is.null(names(x@other)), "elements without names", paste(names(x@other), collapse= "  ")), "\n")
  } else {
    cat("- empty -\n")
  }
  
  cat("\n")
} 
) # end show method for genind




##########################
# Method show for genpop
##########################
setMethod ("show", "genpop", function(object){
  x <- object
  cat("\n")
  cat("       #####################\n")
  cat("       ### Genpop object ### \n")
  cat("       #####################")
  cat("\n- Alleles counts for populations - \n")
  cat("\nS4 class: ", as.character(class(x)))
  
  cat("\n@call: ")
  print(x@call)

  p <- ncol(x@tab)

  cat("\n@tab: ", nrow(x@tab), "x", ncol(x@tab), "matrix of alleles counts\n" )
  
  cat("\n@pop.names: vector of ", length(x@pop.names), "population names")
  cat("\n@loc.names: vector of ", length(x@loc.names), "locus names")
  cat("\n@loc.nall: number of alleles per locus")
  cat("\n@loc.fac: locus factor for the ", ncol(x@tab), "columns of @tab")
  cat("\n@all.names: list of ", length(x@all.names), "components yielding allele names for each locus")
 
  cat("\n\n@other: ")
  if(!is.null(x@other)){
    cat("a list containing: ")
    cat(ifelse(is.null(names(x@other)), "elements without names", paste(names(x@other), collapse= "  ")), "\n")
  } else {
    cat("- empty -\n")
  }
  
  cat("\n")
  
} 
) # end show method for genpop





############################
# Method summary for genind
############################
setMethod ("summary", "genind", function(object, ...){
  x <- object
  if(!inherits(x,"genind")) stop("To be used with a genind object")
  if(is.null(x@pop)){
    x@pop <- factor(rep(1,nrow(x@tab)))
    x@pop.names <- ""
    names(x@pop.names) <- "P1"
  }

  res <- list()

  res$N <- nrow(x@tab)

  res$pop.eff <- as.numeric(table(x@pop))
  names(res$pop.eff) <- names(x@pop.names)

  res$loc.nall <- x@loc.nall

  temp <- genind2genpop(x,quiet=TRUE)@tab

  res$pop.nall <- apply(temp,1,function(r) sum(r!=0,na.rm=TRUE))

  res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab))

  # auxiliary function to compute observed heterozygosity
  temp <- seploc(x,truenames=FALSE,res.type="matrix")
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

  temp <- genind2genpop(x,pop=rep(1,nrow(x@tab)),quiet=TRUE)
  temp <- makefreq(temp,quiet=TRUE)$tab
  temp.names <- colnames(temp)
  temp <- as.vector(temp)
  names(temp) <- temp.names
  temp <- split(temp,x@loc.fac)
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
}) # end summary.genind





############################
# Method summary for genpop
############################
setMethod ("summary", "genpop", function(object, ...){
  x <- object
  if(!inherits(x,"genpop")) stop("To be used with a genpop object")
  
  res <- list()

  res$npop <- nrow(x@tab)

  res$loc.nall <- x@loc.nall

  res$pop.nall <- apply(x@tab,1,function(r) sum(r>0,na.rm=TRUE))

  res$NA.perc <- 100*sum(is.na(x@tab))/prod(dim(x@tab))

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

} 
)# end summary.genpop



#} # end .initAdegenetClasses()






###############
# Methods "is"
###############
is.genind <- function(x){
  res <- ( is(x, "genind") & validObject(x))
  return(res)
}

is.genpop <- function(x){
  res <- ( is(x, "genpop") & validObject(x))
  return(res)
}






###############################################################
###############################################################
# OTHER FUNCTIONS
###############################################################
###############################################################

#########################
# Function genind2genpop
#########################
genind2genpop <- function(x,pop=NULL,missing=c("NA","0","chi2"),quiet=FALSE){

  if(!is.genind(x)) stop("x is not a valid genind object")
  
  if(is.null(x@pop) && is.null(pop)) stop("pop is not provided either in x or in pop")

  missing <- match.arg(missing)

  if(!quiet) cat("\n Converting data from a genind to a genpop object... \n")
  
  # choose pop argument over x@pop
   if(!is.null(pop)) {
    if(length(pop) != nrow(x@tab)) stop("inconsistent length for factor pop")
    # keep levels in order of appearance
    pop <- as.character(pop)
    pop <- factor(pop, levels=unique(pop))
  } else {
    pop <- x@pop
    # keep levels in order of appearance
    pop <- as.character(pop)
    pop <- factor(pop, levels=unique(pop))
    if(!is.null(x@pop.names)) levels(pop) <- x@pop.names # restore real names
  }

  # make generic pop labels, store real pop names
  # pop.names <- levels(pop) ## no longer used

  # tabcount is a matrix pop x alleles, counting alleles per pop
  # *2 to have alleles count
  f1 <- function(v){
    if(all(is.na(v))) return(NA) else return(sum(v,na.rm=TRUE))
  }

  f2 <- function(v){
    if(all(is.na(v)) || sum(v,na.rm=TRUE)==0) return(NA)
    return(v/(sum(v,na.rm=TRUE)))       
  }

  tabcount <- 2* apply(x@tab,2,function(c) tapply(c,pop,f1))
  # restitute matrix class when only one pop
  if(is.null(dim(tabcount))) {
    lab.col <- names(tabcount)
    tabcount <- matrix(tabcount,nrow=1)
    colnames(tabcount) <- lab.col
  }
##   #meancol <- apply(tabcount,2,function(c) mean(c,na.rm=TRUE)) ## no longer used

##   # NA treatment
##   # Treatment when missing='REPLACE':
##   # if allele 'j' of locus 'k' in pop 'i' is missing, replace the NA by a number 'x' so that
##   # the frequency 'x/s' ('s' being the number of observations in 'k' ) equals the frequency 'f'
##   # computed on the whole data (i.e. considering all pop as one)
##   # Then x must verify:
##   # x/s = f(1-f) => x=f(1-f)s
##   #
##   # - eff.pop is a pop x locus matrix giving the corresponding sum of observations (i.e., 's')
##   # - temp is the same table but duplicated for all alleles
##   # - odd.vec is the vector of 'f(1-f)'
##   # - count.replace is a pop x alleles table yielding appropriate replacement numbers (i.e., 'x')

##   if(!is.na(missing) && any(is.na(tabcount))){
##     if(missing==0) tabcount[is.na(tabcount)] <- 0
##     if(toupper(missing)=="REPLACE") {
##     eff.pop <- t(apply(tabcount,1,function(r) tapply(r,x@loc.fac,sum,na.rm=TRUE)))
##     temp <- t(apply(eff.pop,1,function(r) rep(r,table(x@loc.fac))))

##     freq.allpop <- apply(tabcount,2,sum,na.rm=TRUE)
##     freq.allpop <- unlist(tapply(freq.allpop,x@loc.fac,f2))
##     odd.vec <- freq.allpop/(1-freq.allpop)
  
##     count.replace <- t(apply(temp,1,function(r) r*odd.vec))

##     tabcount[is.na(tabcount)] <- count.replace[is.na(tabcount)]
##     }
##   } # end of NA treatment

  
  ## make final object
  temp <- paste(rep(x@loc.names,x@loc.nall),unlist(x@all.names),sep=".")
  colnames(tabcount) <- temp

  prevcall <- match.call()
  
  res <- genpop(tab=tabcount, prevcall=prevcall)
  res@other <- x@other

  if(missing != "NA"){
      res <- na.replace(res, method=missing, quiet=quiet)
  }

  if(!quiet) cat("\n...done.\n\n")

  return(res)
  
} # end genind2genpop





##################
# Methods old2new
##################
setGeneric("old2new",  function(object) standardGeneric("old2new"))

setMethod("old2new", "genind", function(object){
  x <- object
  res <- new("genind")
  theoLength <- 7
  
  res@tab <- as.matrix(x$tab)
  res@ind.names <- as.character(x$ind.names)
  res@loc.names <- as.character(x$loc.names)
  res@loc.nall <- as.integer(x$loc.nall)
  res@loc.fac <- as.factor(x$loc.fac)
  res@all.names <- as.list(x$all.names)
  if(!is.null(x$pop)) {
      res@pop <- as.factor(x$pop)
      theoLength <- theoLength + 1
  }
  if(!is.null(x$pop.names)) {
      res@pop.names <- as.character(x$pop.names)
      theoLength <- theoLength + 1
  }
  res@call <- match.call()

  if(length(object) > theoLength) warning("optional content else than pop and pop.names was not converted")

  return(res)
})


setMethod("old2new", "genpop", function(object){
  x <- object
  res <- new("genpop")

  res@tab <- as.matrix(x$tab)
  res@pop.names <- as.character(x$pop.names)
  res@loc.names <- as.character(x$loc.names)
  res@loc.nall <- as.integer(x$loc.nall)
  res@loc.fac <- as.factor(x$loc.fac)
  res@all.names <- as.list(x$all.names)
  res@call <- match.call()

  if(length(object)>7) warning("optional content was not converted")

  return(res)
})

