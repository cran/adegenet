###########################
#
# Auxiliary functions for
# adegenet objects
#
# T. Jombart
###########################


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
# Function readExt
###################
.readExt <- function(char){
    temp <- as.character(char)
    temp <- unlist(strsplit(char,"[.]"))
    res <- temp[length(temp)]
    return(res)
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





#######################
# Function adegenetWeb
#######################
adegenetWeb <- function(){
    cat("Opening url \"http://adegenet.r-forge.r-project.org/\" ...\n")
    browseURL("http://adegenet.r-forge.r-project.org/")
}





## ############################
## # Function adegenetTutorial
## ############################
## adegenetTutorial <- function(which=c("basics","spca")){
##     which <- match.arg(which)
##     if(which=="general"){
##         url <- "http://adegenet.r-forge.r-project.org/files/adegenet.pdf"
##         cat("\n")
##         cat("  >> Seeking the general tutorial for adegenet.\n")
##         cat("  >> Opening url \"",url,"\".\n ", sep="")
##         cat("\n")
##         browseURL(url)
##     }
##     if(which=="spca"){
##         url <- "http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf"
##         cat("\n")
##         cat("  >> Seeking the sPCA tutorial for adegenet.\n")
##         cat("  >> Opening url \"",url,"\". \n", sep="")
##         cat("\n")
##         browseURL(url)
##     }
## }





############
# checkType
############
##
## WARNING: this does not work with S4 methods
##
checkType <- function(x){
    if(is.character(x)){
        markType <- x
    } else {
        markType <- x@type
    }

    if(markType=="codom") return() # always ok for codominant markers

    currCall <- as.character(sys.call(sys.parent()))[1]
    currFunction <- sub("[[:space:]]*[(].*","",currCall)
    if(currFunction==".local"){
        warning("Current call not found - stopping check (please report this warning).")
        return()
    }

    ## names of functions which are ok for dominant markers
    PAOk <- c("genind","genpop","genind2genpop","summary","df2genind", "genind2df",
                 "truenames","seppop","na.replace","nLoc","scaleGen","spca","selpop")

    PAWarn <- c("df2genind")

    ## function exists but is experimental
    if(currFunction %in% PAWarn){
        msg <- paste(currFunction,"is implemented but experimental presence/absence markers")
        warning(msg)
        return()
    }

    ## function not implemented
    if(! currFunction %in% PAOk){
        msgError <- paste(currFunction,"is not implemented for presence/absence markers")
        stop(msgError)
    } else return() # else, ok.
} # end checkType






##########
## transp
##########
## AUXIL FUNCTION TO USE TRANSPARENT COLORS
transp <- function(col, alpha=.5){
    res <- apply(col2rgb(col),2, function(c) rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
    return(res)
}



##########
## corner
##########
## AUXIL FUNCTION TO ADD LETTER TO A PLOT
corner <- function(text, posi="topleft",  inset=0.1, ...){
    oxpd <- par("xpd")
    on.exit(par(xpd=oxpd))
    par(xpd=TRUE)
    myUsr <- par("usr")
    xrange <- myUsr[1:2]
    yrange <- myUsr[3:4]
    x.size <- abs(diff(xrange))
    y.size <- abs(diff(yrange))
    inset <- rep(inset, length=2)
    x.inset <- inset[1]
    y.inset <- inset[2]

    if(length(grep("top", posi))==1){
        y <- yrange[2] - y.size*y.inset
    } else {
        y <- yrange[1] + y.size*y.inset
    }

    if(length(grep("right", posi))==1){
        x <- xrange[2] - x.size*x.inset
    } else {
        x <- xrange[1] + x.size*x.inset
    }

    text(x, y, lab=text, ...)
}





###########
## num2col
###########
## translate numeric values into colors of a palette
num2col <- function(x, col.pal=heat.colors, reverse=FALSE,
                    x.min=min(x), x.max=max(x), na.col="green"){
    if(any(is.na(x))) warning("NAs detected in x")
    x[x < x.min] <- x.min
    x[x > x.max] <- x.max
    x <- x-x.min # min=0
    x.max <- x.max-x.min # update x.max
    x <- x/x.max # max=1
    x <- round(x*100)
    x[x<=0] <- 1
    if(!reverse) {
        pal <- col.pal(100)
    } else {
        pal <- rev(col.pal(100))
    }

    res <- pal[x]
    res[is.na(res)] <- na.col
    return(res)
}





###########
## fac2col
###########
## translate a factor into colors of a palette
## colors are randomized based on the provided seed
fac2col <- function(x, col.pal=funky, na.col="grey", seed=1){
    ## get factors and levels
    x <- factor(x)
    lev <- levels(x)
    nlev <- length(lev)

    ## get colors corresponding to levels
    set.seed(seed)
    newseed <- round(runif(1,1,1e9))
    on.exit(set.seed(newseed))
    col <- sample(col.pal(nlev))

    ## get output colors
    res <- rep(na.col, length(x))
    res[!is.na(x)] <- col[as.integer(x[!is.na(x)])]

    ## return
    return(res)
}


## pre-defined palettes ##
## mono color
bluepal <- colorRampPalette(c("lightgrey","blue"))
redpal <- colorRampPalette(c("lightgrey","red"))
greenpal <- colorRampPalette(c("lightgrey","green"))

## bi-color
flame <- colorRampPalette(c("gold","red"))

## tri-color
seasun <- colorRampPalette(c("blue","gold","red"))
lightseasun <- colorRampPalette(c("deepskyblue2","gold","red1"))
deepseasun <- colorRampPalette(c("blue2","gold","red2"))

## psychedelic
funky <- colorRampPalette(c("blue","green3","gold","orange","red","brown4","purple"))

