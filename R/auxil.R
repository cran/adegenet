###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################


##############################
# Method truenames for genind
##############################
setGeneric("truenames", function(x) standardGeneric("truenames"))

setMethod("truenames", signature(x="genind"), function(x){
  
  X <- x@tab
  if(!all(x@ind.names=="")) {rownames(X) <- x@ind.names}

  labcol <- rep(x@loc.names,x@loc.nall)
  labcol <- paste(labcol,unlist(x@all.names),sep=".")
  colnames(X) <- labcol

  if(!is.null(x@pop)){
    pop <- x@pop
    levels(pop) <- x@pop.names
    return(list(tab=X,pop=pop))
  }

  return(X)
}
)





##############################
# Method truenames for genpop
##############################
setMethod("truenames",signature(x="genpop"), function(x){

  X <- x@tab
  if(!all(x@pop.names=="")) {rownames(X) <- x@pop.names}

  labcol <- rep(x@loc.names,x@loc.nall)
  labcol <- paste(labcol,unlist(x@all.names),sep=".")
  colnames(X) <- labcol

  return(X)
})




###########################
# Method seploc for genind
###########################
setGeneric("seploc", function(x, ...) standardGeneric("seploc"))

setMethod("seploc", signature(x="genind"), function(x,truenames=TRUE,res.type=c("genind","matrix")){
  
  if(!is.genind(x)) stop("x is not a valid genind object")
  res.type <- match.arg(res.type)
  if(res.type=="genind") { truenames <- TRUE }
  
  temp <- x@loc.fac
  nloc <- length(levels(temp))
  levels(temp) <- 1:nloc

  kX <- list()
  
  for(i in 1:nloc){
    kX[[i]] <- matrix(x@tab[,temp==i],ncol=x@loc.nall[i])

    if(!truenames){
      rownames(kX[[i]]) <- rownames(x@tab)
      colnames(kX[[i]]) <- paste(names(x@loc.names)[i],names(x@all.names[[i]]),sep=".")
    }else{
      rownames(kX[[i]]) <- x@ind.names
      colnames(kX[[i]]) <- paste(x@loc.names[i],x@all.names[[i]],sep=".")
    }
  }

  if(truenames) {
    names(kX) <- x@loc.names
  } else{
    names(kX) <- names(x@loc.names)
  }

  prevcall <- match.call()
  if(res.type=="genind"){
      kX <- lapply(kX, genind, pop=x@pop, prevcall=prevcall)
      for(i in 1:length(kX)){
          kX[[i]]@other <- x@other
      }
  }
  
  return(kX)  
})



###########################
# Method seploc for genpop
###########################
setMethod("seploc", signature(x="genpop"), function(x,truenames=TRUE,res.type=c("genpop","matrix")){
  
  if(!is.genpop(x)) stop("x is not a valid genpop object")
  res.type <- match.arg(res.type)
  if(res.type=="genpop") { truenames <- TRUE }
 
  temp <- x@loc.fac
  nloc <- length(levels(temp))
  levels(temp) <- 1:nloc

  kX <- list()
  
  for(i in 1:nloc){
    kX[[i]] <- matrix(x@tab[,temp==i],ncol=x@loc.nall[i])

    if(!truenames){
      rownames(kX[[i]]) <- rownames(x@tab)
      colnames(kX[[i]]) <- paste(names(x@loc.names)[i],names(x@all.names[[i]]),sep=".")
    }else{
      rownames(kX[[i]]) <- x@pop.names
      colnames(kX[[i]]) <- paste(x@loc.names[i],x@all.names[[i]],sep=".")
    }
  }

  if(truenames) {
    names(kX) <- x@loc.names
  } else{
    names(kX) <- names(x@loc.names)
  }

  prevcall <- match.call()
  if(res.type=="genpop"){
      kX <- lapply(kX, genpop, prevcall=prevcall)
      for(i in 1:length(kX)){
          kX[[i]]@other <- x@other
      }
  }

  return(kX)  
})




#######################
# Function adegenetWeb
#######################
adegenetWeb <- function(){
  cat("Opening url \"http://adegenet.r-forge.r-project.org/\" ...\n")
  browseURL("http://adegenet.r-forge.r-project.org/")
}





###############
# '$' operator
###############
setMethod("$","genind",function(x,name) {
    return(slot(x,name))
})


setMethod("$<-","genind",function(x,name,value) {
   slot(x,name,check=TRUE) <- value
  return(x)
})


setMethod("$","genpop",function(x,name) {
    return(slot(x,name))
})


setMethod("$<-","genpop",function(x,name,value) {
  slot(x,name,check=TRUE) <- value
  return(x)
})





###############
# '[' operator
###############
## genind
setMethod("[","genind",
          function(x, i, j, ..., treatOther=TRUE, drop=FALSE) {

              if (missing(i)) i <- TRUE
              if (missing(j)) j <- TRUE

              pop <- NULL
              if(is.null(x@pop)) { tab <- truenames(x) }
              if(!is.null(x@pop)) {
                  temp <- truenames(x)
                  tab <- temp$tab
                  pop <- temp$pop
                  pop <- factor(pop[i])
              }
             
              prevcall <- match.call()
              tab <- tab[i, j, ...,drop=FALSE]
              
              res <- genind(tab,pop=pop,prevcall=prevcall)

              ## handle 'other' slot
              nOther <- length(x@other)
              namesOther <- names(x@other)
              counter <- 0
              if(treatOther){
                  f1 <- function(obj,n=nrow(x@tab)){
                      counter <<- counter+1
                      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
                          obj <- obj[i,,drop=FALSE]
                      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
                          obj <- obj[i]
                          if(is.factor(obj)) {obj <- factor(obj)}
                      } else {warning(paste("cannot treat the object",namesOther[counter]))}

                      return(obj)
                  } # end f1

                  res@other <- lapply(x@other, f1) # treat all elements
                  
              } # end treatOther
              
              return(res)
          })


## genpop
setMethod("[","genpop", 
          function(x, i, j, ..., treatOther=TRUE, drop=FALSE) {

              if (missing(i)) i <- TRUE
              if (missing(j)) j <- TRUE

              tab <- truenames(x) 
             
              prevcall <- match.call()
              tab <- tab[i, j, ...,drop=FALSE]
              
              res <- genpop(tab,prevcall=prevcall)

              ## handle 'other' slot
              nOther <- length(x@other)
              namesOther <- names(x@other)
              counter <- 0
              if(treatOther){
                  f1 <- function(obj,n=nrow(x@tab)){
                      counter <<- counter+1
                      if(!is.null(dim(obj)) && nrow(obj)==n) { # if the element is a matrix-like obj
                          obj <- obj[i,,drop=FALSE]
                      } else if(length(obj) == n) { # if the element is not a matrix but has a length == n
                          obj <- obj[i]
                          if(is.factor(obj)) {obj <- factor(obj)}
                      } else {warning(paste("cannot treat the object",namesOther[counter]))}
                      
                      return(obj)
                  } # end f1
                  
                  res@other <- lapply(x@other, f1) # treat all elements
                  
              } # end treatOther
             
              
              return(res)
          })






##################
# Function seppop
##################
setGeneric("seppop", function(x, ...) standardGeneric("seppop"))

## genind
setMethod("seppop", signature(x="genind"), function(x,pop=NULL,truenames=TRUE,res.type=c("genind","matrix")){

    ## misc checks 
    if(!is.genind(x)) stop("x is not a valid genind object")
    if(is.null(pop)) {pop <- x@pop}
    if(is.null(pop)) stop("pop not provided and x@pop is empty")
    res.type <- match.arg(res.type)
    if(res.type=="genind") { truenames <- TRUE }
  
    pop <- x@pop
    levels(pop) <- x@pop.names

    ## make a list of genind objects
    kObj <- lapply(levels(pop), function(lev) x[pop==lev, ])
    names(kObj) <- levels(pop)

    ## res is a list of genind
    if(res.type=="genind"){ return(kObj) }
  
    ## res is list of matrices
    if(truenames) {
        res <- lapply(kObj, function(obj) truenames(obj)$tab)
    } else{
        res <- lapply(kObj, function(obj) obj$tab)
    }
    
    return(res) 
}) # end seppop





#####################
# Methods na.replace
#####################
setGeneric("na.replace", function(x, ...) standardGeneric("na.replace"))

## genind method
setMethod("na.replace", signature(x="genind"), function(x,method, quiet=FALSE){

    ## preliminary stuff
    validObject(x)
    if(!any(is.na(x@tab))) {
        if(!quiet) cat("\n Replaced 0 missing values \n")
        return(x)
    }
    method <- tolower(method)
    method <- match.arg(method, c("0","mean"))

    res <- x
    
    if(method=="0"){
        res@tab[is.na(x@tab)] <- 0
    }

    if(method=="mean"){
        f1 <- function(vec){
            m <- mean(vec,na.rm=TRUE)
            vec[is.na(vec)] <- m
            return(vec)
        }

        res@tab <- apply(x@tab, 2, f1)
    }

    if(!quiet){
        Nna <- sum(is.na(x@tab))
        cat("\n Replaced",Nna,"missing values \n")
    }

    return(res)

})




## genpop method
setMethod("na.replace", signature(x="genpop"), function(x,method, quiet=FALSE){

    ## preliminary stuff
    validObject(x)
    if(!any(is.na(x@tab))) {
        if(!quiet) cat("\n Replaced 0 missing values \n")
        return(x)
    }

    method <- tolower(method)
    method <- match.arg(method, c("0","chi2"))

    res <- x
    
    if(method=="0"){
        res@tab[is.na(x@tab)] <- 0
    }

    if(method=="chi2"){
        ## compute theoretical counts
        ## (same as in a Chi-squared)
        X <- x@tab
        sumPop <- apply(X,1,sum,na.rm=TRUE)
        sumLoc <- apply(X,2,sum,na.rm=TRUE)
        X.theo <- sumPop %o% sumLoc / sum(X,na.rm=TRUE)

        X[is.na(X)] <- X.theo[is.na(X)]
        res@tab <- X
    }

    if(!quiet){
        Nna <- sum(is.na(x@tab))
        cat("\n Replaced",Nna,"missing values \n")
    }

    return(res)
})





##################
# Function repool
##################
repool <- function(...){

    ## preliminary stuff
    x <- list(...)
    if(is.list(x[[1]])) x <- x[[1]] ## if ... is a list, keep this list for x
    if(!inherits(x,"list")) stop("x must be a list")
    if(!all(sapply(x,is.genind))) stop("x is does not contain only valid genind objects")
    temp <- sapply(x,function(e) e$loc.names)
    if(!all(table(temp)==length(x))) stop("markers are not the same for all objects")
    
    
    ## extract info
    listTab <- lapply(x,genind2df,usepop=FALSE)
    getPop <- function(obj){
        if(is.null(obj$pop)) return(factor(rep(NA,nrow(obj$tab))))
      pop <- obj$pop
        levels(pop) <- obj$pop.names
        return(pop)
    }
    
    ## handle pop
    listPop <- lapply(x, getPop)
    pop <- unlist(listPop, use.name=FALSE)
    pop <- factor(pop)
    
  ## handle genotypes
    markNames <- colnames(listTab[[1]])
    listTab <- lapply(listTab, function(tab) tab[,markNames]) # resorting of the tabs
    
    ## bind all tabs by rows
    tab <- listTab[[1]] 
    for(i in 2:length(x)){
        tab <- rbind(tab,listTab[[i]])
    }
    
    res <- df2genind(tab,pop=pop)
    res$call <- match.call()
    
    return(res)
} # end repool


