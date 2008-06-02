## propShared computes the proportion of shared alleles
## in a genind object


######################
# Function propShared
######################
propShared <- function(obj){
    x <- obj
    ## check that this is a valid genind
    if(!inherits(x,"genind")) stop("obj must be a genind object.")
    invisible(validObject(x))
    
    ## build a matrix of genotypes (in rows) coded by integers
    ## NAs are coded by 0
    ## The matrix is a cbind of two matrices, storing respectively the
    ## first and the second allele.
    temp <- genind2df(x,usepop=FALSE)
    alleleSize <- max(apply(temp,1:2,nchar))/2
    mat1 <- apply(temp, 1:2, substr, 1, alleleSize)
    mat2 <- apply(temp, 1:2, substr, alleleSize+1, alleleSize*2)
    matAll <- cbind(mat1,mat2)
    matAll <- apply(matAll,1:2,as.integer)
    matAll[is.na(matAll)] <- 0

    n <- nrow(matAll)
    resVec <- double(n*(n-1)/2)
    res <- .C("sharedAll", as.integer(as.matrix(matAll)),
              n, ncol(matAll), resVec, PACKAGE="adegenet")[[4]]

    attr(res,"Size") <- n
    attr(res,"Diag") <- FALSE
    attr(res,"Upper") <- FALSE
    class(res) <- "dist"
    res <- as.matrix(res)

    diag(res) <- 1
    rownames(res) <- x@ind.names
    colnames(res) <- x@ind.names
    return(res)
}




## ######################
## # Function propShared (old, pure-R version)
## ######################
## propShared <- function(obj){

##     x <- obj
##     ## check that this is a valid genind
##     if(!inherits(x,"genind")) stop("obj must be a genind object.")
##     invisible(validObject(x))
    
##     ## replace NAs
##     x <- na.replace(x, method="0")
    
##     ## some useful variables
##     nloc <- length(x@loc.names)
    
##     ## fnorm: auxiliary function for scaling
##     fnorm <- function(vec){
## 	norm <- sqrt(sum(vec*vec))
## 	if(length(norm) > 0 && norm > 0) {return(vec/norm)}
## 	return(vec)	
##     }
    
##     ## auxiliary function f1
##     ## computes the proportion of shared alleles in one locus
##     f1 <- function(X){
## 	X <- t(apply(X, 1, fnorm))
## 	res <- X %*% t(X)
## 	res[res>0.51 & res<0.9] <- 0.5 # remap case one heteroZ shares the allele of on homoZ
## 	return(res)
##     }
    
##     ## separate data per locus
##     temp <- seploc(x)
##     listProp <- lapply(temp, function(e) f1(e@tab))
    
##     ## produce the final result
##     res <- listProp[[1]]
##     if(nloc>2){
## 	for(i in 2:nloc){
##             res <- res + listProp[[i]]
## 	}
##     }
    
##     res <- res/nloc
##     rownames(res) <- x@ind.names
##     colnames(res) <- x@ind.names
    
##     return(res)
## }
