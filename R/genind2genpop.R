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
  # *ploidy to have alleles counts
  f1 <- function(v){
    if(all(is.na(v))) return(NA) else return(sum(v,na.rm=TRUE))
  }

  f2 <- function(v){
    if(all(is.na(v)) || sum(v,na.rm=TRUE)==0) return(NA)
    return(v/(sum(v,na.rm=TRUE)))
  }

  tabcount <- x@ploidy * apply(x@tab,2,function(c) tapply(c,pop,f1))
  tabcount <- round(tabcount,digits=0)
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
