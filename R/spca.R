##############################################
#
# spatial Principal Components Analysis
#
# require ade4, spdep and eventually tripack
#
# generic functions were derived from
# those of multispati class (ade4)
# 
# T. Jombart (jombart@biomserv.univ-lyon1.fr)
# 31 may 2007
##############################################



################
# Function spca
################
spca <- function(obj, xy=NULL, cn=NULL, scale=FALSE, scannf=TRUE, nfposi=1, nfnega=1, type=NULL,
                 ask=TRUE, plot.nb=TRUE, edit.nb=FALSE ,truenames=TRUE, d1=NULL, d2=NULL, k=NULL,
                 a=NULL, dmin=NULL) {
  
  if(!any(inherits(obj,c("genind","genpop")))) stop("obj must be a genind or genpop object.")
  invisible(validObject(obj))

  ## spatial coordinates
  if(is.null(xy) & !is.null(obj$other$xy)) xy <- obj$other$xy
  if(is.data.frame(xy)) xy <- as.matrix(xy)
  if(!is.null(xy) & !is.matrix(xy)) stop("wrong 'xy' provided")
  
  if(!require(ade4, quiet=TRUE)) stop("ade4 library is required.")

  appel <- match.call()
  
  ## connection network
  if(is.null(cn)) {
    if(is.null(xy)) stop("'xy' and 'cn' are both missing")
    resCN <- chooseCN(xy=xy, ask=ask, type=type, plot.nb=plot.nb, edit.nb=edit.nb,
                      result.type="listw", d1=d1, d2=d2, k=k, a=a, dmin=dmin)
  } else {
      if(inherits(cn,"nb") & !inherits(cn,"listw")) {
          xy <- attr(cn,"xy") # xy coords can be retrieved from cn of class nb (not from listw) 
          cn <- nb2listw(cn, style="W", zero.policy=TRUE)
      }

      if(!inherits(cn,"listw")) {
          stop("cn does not have a recognized class ('nb' or 'listw', package spdep)")
      } else {
          if(is.null(xy)) stop("listw object provided as 'cn' without providing 'xy'")
          resCN <- cn
      }
  }

  ## check xy coordinates
  if(ncol(xy) != 2) stop("xy does not have two columns.")
  if(nrow(xy) != nrow(obj@tab)) stop("obj@tab and xy must have the same row numbers.")

  ## prepare data
  f1 <- function(vec){
    m <- mean(vec,na.rm=TRUE)
    vec[is.na(vec)] <- m
    return(vec)
  }

  if(is.genind(obj)) { X <- obj@tab }
  if(is.genpop(obj)) { X <- makefreq(obj, quiet=TRUE)$tab }

  ## handle NAs
  if(any(is.na(X))){
      warning("NAs in data are automatically replaced (to mean allele frequency")
      X <- apply(X,2,f1)
  }

  if(truenames){
    rownames(X) <- rownames(truenames(obj))
    colnames(X) <- colnames(truenames(obj))   
  }

  # perform analyses
  pcaX <- dudi.pca(X, center=TRUE, scale=scale, scannf=FALSE)

  spcaX <- multispati(dudi=pcaX, listw=resCN, scannf=scannf, nfposi=nfposi, nfnega=nfnega)

  nfposi <- spcaX$nfposi
  nfnega <- spcaX$nfnega

  spcaX$xy <- xy
  rownames(spcaX$xy) <- rownames(spcaX$li)
  colnames(spcaX$xy) <- c("x","y")
  
  spcaX$lw <- resCN
  
  spcaX$call <- appel

  posaxes <- if(nfposi>0) {1:nfposi} else NULL
  negaxes <- if(nfnega>0) {(length(spcaX$eig)-nfnega+1):length(spcaX$eig)} else NULL
  keptaxes <- c(posaxes,negaxes)
      
  colnames(spcaX$c1) <- paste("Axis",keptaxes)
  colnames(spcaX$li) <- paste("Axis",keptaxes)
  colnames(spcaX$ls) <- paste("Axis",keptaxes)

  class(spcaX) <- "spca"

  return(spcaX)

} # end spca





######################
# Function print.spca
######################
print.spca <- function(x, ...){
  cat("\t########################################\n")
  cat("\t# Spatial principal component analysis #\n")
  cat("\t########################################\n")
  cat("class: ")
  cat(class(x))
  cat("\n$call: ")
  print(x$call)
  cat("\n$nfposi:", x$nfposi, "axis-components saved")
  cat("\n$nfnega:", x$nfnega, "axis-components saved")
 
  cat("\nPositive eigenvalues: ")
  l0 <- sum(x$eig >= 0)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5) 
    cat(" ...\n")
  else cat("\n")  
  cat("Negative eigenvalues: ")
  l0 <- sum(x$eig <= 0)
  cat(sort(signif(x$eig, 4))[1:(min(5, l0))])
  if (l0 > 5) 
    cat(" ...\n")
  else cat("\n")
  cat('\n')
  sumry <- array("", c(1, 4), list(1, c("vector", "length", 
                                        "mode", "content")))
  sumry[1, ] <- c('$eig', length(x$eig), mode(x$eig), 'eigenvalues')
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(4, 4), list(1:4, c("data.frame", "nrow", "ncol", "content")))
  sumry[1, ] <- c("$c1", nrow(x$c1), ncol(x$c1), "principal axes: scaled vectors of alleles loadings")
  sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "principal components: coordinates of entities")
  sumry[3, ] <- c("$ls", nrow(x$ls), ncol(x$ls), 'lag vector of principal components')
  sumry[4, ] <- c("$as", nrow(x$as), ncol(x$as), 'pca axes onto spca axes')
  
  class(sumry) <- "table"
  print(sumry)

  cat("\n$xy: matrix of spatial coordinates")
  cat("\n$lw: a list of spatial weights (class 'listw')")
  
  cat("\n\nother elements: ")
  if (length(names(x)) > 10) 
    cat(names(x)[11:(length(names(x)))], "\n")
  else cat("NULL\n")
}





########################
# Function summary.spca
########################
summary.spca <- function (object, ..., printres=TRUE) {
  if (!inherits(object, "spca"))stop("to be used with 'spca' object")
  if(!require(ade4,quietly=TRUE)) stop("the library ade4 is required; please install this package")
  if(!require(spdep,quietly=TRUE)) stop("the library spdep is required; please install this package")

  #util <- function(n) { ## no longer used
  #  x <- "1"
  #  for (i in 2:n) x[i] <- paste(x[i - 1], i, sep = "+")
  #  return(x)
  #}
  norm.w <- function(X, w) {
    f2 <- function(v) sum(v * v * w)/sum(w)
    norm <- apply(X, 2, f2)
    return(norm)
  }

  resfin <- list()
  
  if(printres) {
    cat("\nSpatial principal component analysis\n")
    cat("\nCall: ")
    print(object$call)
  }
  
  appel <- as.list(object$call)
  ## compute original pca
  # prepare data
  obj <- eval(appel$obj)
  if(is.null(appel$truenames)) truenames <- FALSE
  
  f1 <- function(vec){
    m <- mean(vec,na.rm=TRUE)
    vec[is.na(vec)] <- m
    return(vec)
  }
  
  if(is.genind(obj)) { X <- obj@tab }
  if(is.genpop(obj)) { X <- makefreq(obj, quiet=TRUE)$tab }
  
  X <- apply(X,2,f1)
  
  if(truenames){
    rownames(X) <- rownames(truenames(obj))
    colnames(X) <- colnames(truenames(obj))   
  }
  
  nfposi <- object$nfposi
  nfnega <- object$nfnega
  
  dudi <- dudi.pca(X, center=TRUE, scale=FALSE, scannf=FALSE, nf=nfposi+nfnega)
  ## end of pca
    
  lw <- object$lw

  # I0, Imin, Imax
  n <- nrow(X)
  I0 <- -1/(n-1)
  L <- listw2mat(lw)
  eigL <- eigen(0.5*(L+t(L)))$values
  Imin <- min(eigL)
  Imax <- max(eigL)
  Ival <- data.frame(I0=I0,Imin=Imin,Imax=Imax)
  row.names(Ival) <- ""
  if(printres) {
    cat("\nConnection network statistics:\n")
    print(Ival)
  }

  Istat <- c(I0,Imin,Imax)
  names(Istat) <- c("I0","Imin","Imax")
  resfin$Istat <- Istat

  
  # les scores de l'analyse de base
  nf <- dudi$nf
  eig <- dudi$eig[1:nf]
  cum <- cumsum(dudi$eig)[1:nf]
  ratio <- cum/sum(dudi$eig)
  w <- apply(dudi$l1,2,lag.listw,x=lw)
  moran <- apply(w*as.matrix(dudi$l1)*dudi$lw,2,sum)
  res <- data.frame(var=eig,cum=cum,ratio=ratio, moran=moran)
  row.names(res) <- paste("Axis",1:nf)
  if(printres) {
    cat("\nScores from the centred PCA\n")
    print(res)
  }

  resfin$pca <- res

  
  # les scores de l'analyse spatiale
  # on recalcule l'objet en gardant tous les axes
  eig <- object$eig
  nfposimax <- sum(eig > 0)
  nfnegamax <- sum(eig < 0)
    
  ms <- multispati(dudi=dudi, listw=lw, scannf=FALSE,
                     nfposi=nfposimax, nfnega=nfnegamax)

  ndim <- dudi$rank
  nf <- nfposi + nfnega
  agarder <- c(1:nfposi,if (nfnega>0) (ndim-nfnega+1):ndim)
  varspa <- norm.w(ms$li,dudi$lw)
  moran <- apply(as.matrix(ms$li)*as.matrix(ms$ls)*dudi$lw,2,sum)
  res <- data.frame(eig=eig,var=varspa,moran=moran/varspa)
  row.names(res) <- paste("Axis",1:length(eig))
  
  if(printres) {
    cat("\nsPCA eigenvalues decomposition:\n")
    print(res[agarder,])
  }
  
  resfin$spca <- res
    
  return(invisible(resfin))
}



#####################
# Function plot.spca
#####################
plot.spca <- function (x, axis = 1, ...){
    if (!inherits(x, "spca")) stop("Use only with 'spca' objects.")

    if(!require(ade4)) stop("ade4 package is required.")
    if(!require(spdep)) stop("spdep package is required.")
    if(axis>ncol(x$li)) stop("wrong axis required.")
    
    opar <- par(no.readonly = TRUE)
    on.exit(par(opar))
    par(mar = rep(.1,4), mfrow=c(3,2))

    n <- nrow(x$li)
    xy <- x$xy
    z <- x$ls[,axis]
    nfposi <- x$nfposi
    nfnega <- x$nfnega
    ## handle neig parameter - hide cn if nore than 100 links
    nLinks <- sum(card(x$lw$neighbours))
    if(nLinks < 500) {
        neig <- nb2neig(x$lw$neighbours)
    } else {
        neig <- NULL
    }
    
    sub <- paste("Score",axis)
    csub <- 2
      
    # 1
    if(n<30) clab <- 1 else clab <- 0
    s.label(xy, clab=clab, include.ori=FALSE, addaxes=FALSE, neig=neig,
            cneig=1, sub="Connection network", csub=2)    
    
    # 2
    s.image(xy,z, include.ori=FALSE, grid=TRUE, kgrid=10, cgrid=1,
            sub=sub, csub=csub, possub="bottomleft")
    box()
    
    # 3
    if(n<30) {neig <- nb2neig(x$lw$neighbours)} else {neig <- NULL}
    s.value(xy,z, include.ori=FALSE, addaxes=FALSE, clegend=0, csize=.6,
            neig=neig, sub=sub, csub=csub, possub="bottomleft")
    
    # 4
    s.value(xy,z, include.ori=FALSE, addaxes=FALSE, clegend=0, csize=.6,
            method="greylevel", neig=neig, sub=sub, csub=csub, possub="bottomleft")
        
    # 5
    omar <- par("mar")
    par(mar = c(0.8, 2.8, 0.8, 0.8))
    m <- length(x$eig)
    col.w <- rep("white", m) # elles sont toutes blanches
    col.w[1:nfposi] <- "grey"
    if (nfnega>0) {col.w[m:(m-nfnega+1)] <- "grey"}
    j <- axis
    if (j>nfposi) {j <- j-nfposi +m -nfnega}
    col.w[j] <- "black" 
    barplot(x$eig, col = col.w)
    scatterutil.sub(cha ="Eigenvalues", csub = 2.5, possub = "topright")
    par(mar=rep(.1,4))
    box()
    par(mar=omar)
    
    # 6
    par(mar=c(4,4,2,1))
    screeplot(x,main="Eigenvalues decomposition")
    par(mar=rep(.1,4))
    box()
    return(invisible(match.call()))
}



##########################
# Function screeplot.spca
##########################
screeplot.spca <- function(x,...,main=NULL){

  opar <- par("las")
  on.exit(par(las=opar))

  sumry <- summary(x,printres=FALSE)
  
  labels <- lapply(1:length(x$eig),function(i) bquote(lambda[.(i)]))

  par(las=1)
  
  xmax <- sumry$pca[1,1]*1.1
  I0 <- sumry$Istat[1]
  Imin <- sumry$Istat[2]
  Imax <- sumry$Istat[3]

  plot(x=sumry$spca[,2],y=sumry$spca[,3],type='n',xlab='Variance',ylab="Spatial autocorrelation (I)",xlim=c(0,xmax),ylim=c(Imin*1.1,Imax*1.1),yaxt='n',...)
  text(x=sumry$spca[,2],y=sumry$spca[,3],do.call(expression,labels))
  
  ytick <- c(I0,round(seq(Imin,Imax,le=5),1))
  ytlab <- as.character(round(seq(Imin,Imax,le=5),1))
  ytlab <- c(as.character(round(I0,1)),as.character(round(Imin,1)),ytlab[2:4],as.character(round(Imax,1)))
  axis(side=2,at=ytick,labels=ytlab)

  rect(0,Imin,xmax,Imax,lty=2)
  segments(0,I0,xmax,I0,lty=2)
  abline(v=0)

  if(is.null(main)) main <- ("Spatial and variance components of the eigenvalues")
  title(main)
  
  return(invisible(match.call()))
}
