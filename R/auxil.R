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





############################
# Function adegenetTutorial
############################
adegenetTutorial <- function(which=c("general","spca")){
    which <- match.arg(which)
    if(which=="general"){
        url <- "http://adegenet.r-forge.r-project.org/files/adegenet.pdf"
        cat("\n")
        cat("  >> Seeking the general tutorial for adegenet.\n")
        cat("  >> Opening url \"",url,"\".\n ", sep="")
        cat("\n")
        browseURL(url)
    }
    if(which=="spca"){
        url <- "http://adegenet.r-forge.r-project.org/files/tutorial-spca.pdf"
        cat("\n")
        cat("  >> Seeking the sPCA tutorial for adegenet.\n")
        cat("  >> Opening url \"",url,"\". \n", sep="")
        cat("\n")
        browseURL(url)
    }
}





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


