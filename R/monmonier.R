# Algorithm to detect boundaries, based on Monmonier's algorithm
# Extended to any connection network 
# Thibaut Jombart 2006-2007 (jombart@biomserv.univ-lyon1.fr)



#####################
# function monmonier
#####################
monmonier <- function(xy,dist,cn,threshold=NULL,nrun=1,skip.local.diff=rep(0,nrun),scanthres=is.null(threshold)){
if(!require(spdep) & !require(ade4)) stop("The package spdep is required but not installed")
if(!inherits(cn,"nb")) stop('cn is not a nb object')
if(is.data.frame(xy)) xy <- as.matrix(xy)
if(!is.matrix(xy)) stop('xy must be a matrix')
if(!inherits(dist,"dist")) stop('Argument \'dist\' must be a distance matrix of class dist')
if(nrow(xy) != nrow(as.matrix(dist))) stop('Number of sites and number of observations differ')

# precision des coordonnees xy (nb chiffres apres virgule)
# ne pas dépasser 10 !
PRECISION=8

# conversion cn
cn.nb <- cn
cn <- nb2neig(cn)
# calcul matrice voisinage
M <- neig2mat(cn)
# matrice de distance D
D <- as.matrix(dist)
# matrice des distances entre voisins, D
D <- M*D
# mettre la valeur seuil, par expl la médiane des différences entre voisins.
if(is.null(threshold) || threshold<0) {Dlim <- summary(unique(D[D>0]))[5]} else {Dlim <- threshold}

if(scanthres){
  barplot(sort(unique(D)[unique(D) > 0],dec=TRUE),main="Local distances barplot")
  abline(h=Dlim,lty=2)
  mtext("Dashed line indicates present threshold")
  cat("Indicate the threshold (\'d\' for default): ")
  temp <- as.character(readLines(n = 1))
  if(toupper(temp)!="D") { Dlim <- as.numeric(temp) }
}

# faire liste des coord des points voisins, data.frame avec comme colonnes : x1 y1 x2 y2.
listCpl <- neig.util.GtoL(neig2mat(cn))
allSeg <- cbind(xy[,1][listCpl[,1]] , xy[,2][listCpl[,1]] , xy[,1][listCpl[,2]] , xy[,2][listCpl[,2]])
colnames(allSeg) <- c('xP','yP','xQ','yQ')

# retourne la valeur d'une différence selon le rang (ordre décroissant) et les numéros des points correspondants
# et retourne les coord du milieu M d'un segment AB.
getNext <- function(D,rang){
	val <- unique(sort(D,decreasing=TRUE))[rang]
	A <- which(D==val,TRUE)[1,1]
	B <- which(D==val,TRUE)[1,2]
	xA <- xy[A,1]
	yA <- xy[A,2]
	xB <- xy[B,1]
	yB <- xy[B,2]
	xM <- (xA + xB)/2
	yM <- (yA + yB)/2
	return( list(A=c(A,xA,yA), B=c(B,xB,yB), M=c(xM,yM), val=val) )
}

## retourne FALSE si on croise avec une arrete, sinon TRUE
# M et N definissent le segment d'interet
# segMat est la matrice des segments connus (arretes)
# AB est le segment de milieu M ; on le donne pour enlever la comparaison MN vs AB
# CD est le segment de milieu N ; on le donne pour enlever la comparaison MN vs CD
checkNext <- function(M,N,A,B,C,D,segMat=allSeg){
        # orientation du segment MN
	xmin <- min(M[1],N[1])
	xmax <- max(M[1],N[1])
	ymin <- min(M[2],N[2])
	ymax <- max(M[2],N[2])

        # A partir d'ici on élimine des comparaisons inutiles et/ou sources de problèmes
        # Faire attention à toujours garder subsetSeg en tant que matrice
        
        # il faut éliminer le segment qui vient d'être tracé, des fois que les comparaisons XY vs XZ renvoient un code d'intersection ordinaire.
        subsetSeg <- segMat[-nrow(segMat),]
        subsetSeg <- matrix(subsetSeg,ncol=4) # coerce to matrix, even if 0 rows
        
        # subsetSeg est une matrice dont chaque ligne est un segment : xP,yP,xQ,yQ
        # on elimine les segments totalement hors du carre de diagonale MN
	subsetSeg <- matrix(subsetSeg[!(subsetSeg[,1] < xmin & subsetSeg[,3] < xmin),] ,ncol=4)
	subsetSeg <- matrix(subsetSeg[!(subsetSeg[,1] > xmax & subsetSeg[,3] > xmax),] ,ncol=4)
	subsetSeg <- matrix(subsetSeg[!(subsetSeg[,2] < ymin & subsetSeg[,4] < ymin),] ,ncol=4)
	subsetSeg <- matrix(subsetSeg[!(subsetSeg[,2] > ymax & subsetSeg[,4] > ymax),] ,ncol=4)
        # ntemp est le nombre de ligne de subsetSeg a ce stade
        # a partir d'ici, on va eliminer AB et CD ; donc si ntemp <= 2, pas la peine de continuer
        ntemp <- nrow(subsetSeg)
        if(ntemp <= 2) return(TRUE)
        
        # on elimine le segment AB dont le milieu est N (code 2 litigieux)
        # il faut le retrouver dans la liste des segments
        # on compare pour se faire les coordonnees, arrondies puisqu'en double
        # il faut aussi savoir ou sont A et B (en premier ou deuxieme) dans le tableau
        AB <- c(A,B)
        temp <- apply(subsetSeg,1,function(r) all(round(r,PRECISION)==round(AB,PRECISION)) )
        if(!any(temp)) {
          AB <- c(B,A)
          temp <- apply(subsetSeg,1,function(r) all(round(r,PRECISION)==round(AB,PRECISION)) )
        }
        if(!any(temp)) {
          # This warning is no longer useful. Commented.   
          # warning("Failed to avoid middle-segment comparaison. Wrong path is likely to result.")
        } else{ subsetSeg <- subsetSeg[-which(temp),] }

        # idem pour CD, qu'il faut enlever      
        CD <- c(C,D)
        temp <- apply(subsetSeg,1,function(r) all(round(r,PRECISION)==round(CD,PRECISION)) )
        if(!any(temp)) {
          CD <- c(D,C)
          temp <- apply(subsetSeg,1,function(r) all(round(r,PRECISION)==round(CD,PRECISION)) )
        }
        if(!any(temp)) {
          # This warning is no longer useful. Commented.    
          # warning("Failed to avoid middle-segment comparaison. Wrong path is likely to result.")}
        } else{ subsetSeg <- subsetSeg[-which(temp),] }

        # temp utilisé pour la réponse (ecrase l'ancien), initialisé à 10
	temp <- as.integer(10)
                
        # version utilisant monmonier-utils.C
        # evalue seulement si il y a des comparaisons a faire
        # attention : si ntemp=3 (avant de virer AB et CD), alors ici subsetSeg n'est plus une matrice, mais un vecteur
        if(ntemp == 3) {subsetSeg <- matrix(subsetSeg,ncol=4)}
        if(nrow(subsetSeg)>0) {
          temp <- .C("CheckAllSeg",as.integer(nrow(subsetSeg)),as.integer(ncol(subsetSeg)),
          as.double(as.matrix(subsetSeg)), as.double(M), as.double(N), temp,PACKAGE="adegenet")[[6]]
        } else {temp <- 0}
        
	# chargement de la fonction C (a commenter une fois dans ade4)
	#dyn.load('monmonier-utils.so') !!! marche pas, besoin de taballoc

	# si 1 (au moins une intersection) est retourné, on retourne FALSE, TRUE dans tous les autres cas.
        # on ajoute un controle
        if(temp==10) stop("CheckAllSeg failure (returned value=10, i.e. unchanged, not computed). Please report the bug.")
        if(temp==1) return(FALSE) else return(TRUE)
}


# création de l'objet result
# c'est une liste ayant une liste de résultats par run
result <-list()
for(run in 1:nrun){
	result[[run]] <- list(dir1=list(),dir2=list())
}

# initialiser en trouvant la plus grande différence entre voisin
# puis tracer le chemin en fonction du premier point
# se fait en prenant la plus grande différence pour le run 1, la seconde pour le 2, etc.
for(run in 1:nrun){
  
        ## dir 1 ##

        # on retire éventuellement les différences locales qui piegent l'algorithme (argument skip)
	# ces valeurs sont retirées de facon irrémédiables (mises à -1)
	if(skip.local.diff[run] >0) 
		for(i in 1:skip.local.diff[run]){
			temp <- getNext(D,1)
			D[temp$A[1],temp$B[1] ] <- -1
			D[temp$B[1],temp$A[1] ] <- -1
		}
	# initialisation
	temp <- getNext(D,1)
        firstPoint <- temp
	if(temp$val<=Dlim) stop(paste('Algorithm reached the threshold value at the first step of run',run))
	result[[run]]$dir1[[1]] <- temp
	D[result[[run]]$dir1[[1]]$A[1],result[[run]]$dir1[[1]]$B[1]] <- -1
	D[result[[run]]$dir1[[1]]$B[1],result[[run]]$dir1[[1]]$A[1]] <- -1

	### fonction principale
	# i est l'indice du résultat, s est l'indice de rang de la différence entre voisin (ordre decroissant)
	i <- 1
	s <- 1
	while(temp$val>Dlim){
          temp <- getNext(D,s) 
		if( (checkNext(result[[run]]$dir1[[length(result[[run]]$dir1)]]$M,
                               temp$M,
                               result[[run]]$dir1[[length(result[[run]]$dir1)]]$A[2:3],
                               result[[run]]$dir1[[length(result[[run]]$dir1)]]$B[2:3],
                               temp$A[2:3],
                               temp$B[2:3])) & (temp$val>Dlim)){
			i <- i+1
			result[[run]]$dir1[[i]] <- temp
                        # mise à jour matrice des différences
			D[result[[run]]$dir1[[i]]$A[1],result[[run]]$dir1[[i]]$B[1]] <- -1
			D[result[[run]]$dir1[[i]]$B[1],result[[run]]$dir1[[i]]$A[1]] <- -1
			s <- 1
			# mise a jour des segments
			allSeg <- rbind(allSeg,c(result[[run]]$dir1[[i-1]]$M,result[[run]]$dir1[[i]]$M) ) 
		
		} else{ s <- s+1 }
        
	} # end dir 1 pour un run donne
        
        ## dir 2 ##
        temp <- firstPoint
          
        result[[run]]$dir2[[1]] <- temp
        D[result[[run]]$dir2[[1]]$A[1],result[[run]]$dir2[[1]]$B[1]] <- -1
        D[result[[run]]$dir2[[1]]$B[1],result[[run]]$dir2[[1]]$A[1]] <- -1
        
        # i est l'indice du résultat, s est l'indice de rang de la différence entre voisin (ordre decroissant)
        i <- 1
        s <- 1
        while(temp$val>Dlim){
	  temp <- getNext(D,s)
          if( checkNext(result[[run]]$dir2[[length(result[[run]]$dir2)]]$M,
                        temp$M,
                        result[[run]]$dir2[[length(result[[run]]$dir2)]]$A[2:3],
                        result[[run]]$dir2[[length(result[[run]]$dir2)]]$B[2:3],
                        temp$A[2:3],
                        temp$B[2:3]) & (temp$val>Dlim)){
		   i <- i+1
                   result[[run]]$dir2[[i]] <- temp
                   D[result[[run]]$dir2[[i]]$A[1],result[[run]]$dir2[[i]]$B[1]] <- -1
                   D[result[[run]]$dir2[[i]]$B[1],result[[run]]$dir2[[i]]$A[1]] <- -1
                   s <- 1
                   # mise a jour des segments
                   allSeg <- rbind(allSeg,c(result[[run]]$dir2[[i-1]]$M,result[[run]]$dir2[[i]]$M) )		
          } else{ s <- s+1 }        
              } # end dir2 for a given run
} # end for all run

# mise en forme de l'output
# c'est une liste de la classe monmonier
# chaque élément correspond à un run, donc à une "frontière" potentielle.
# l'objet contient aussi le nombre de runs ($nrun) et son propre appel ($call).
output=list()
for(run in 1:nrun){
	runname <- paste('run',run,sep='')
	output[[runname]] <- list(dir1=list(),dir2=list())
	# dir 1 #
	output[[runname]]$dir1$path <- matrix(-1, ncol=2,nrow=length(result[[run]]$dir1))
	colnames(output[[runname]]$dir1$path) <- c('x','y')
	rownames(output[[runname]]$dir1$path) <- paste('Point',1:nrow(output[[run]]$dir1$path),sep='_')

	for(i in 1:length(result[[run]]$dir1)) {
		output[[runname]]$dir1$path[i,] <- result[[run]]$dir1[[i]]$M
		output[[runname]]$dir1$values[i] <- result[[run]]$dir1[[i]]$val
	}

        # dir 2 #
	output[[runname]]$dir2$path <- matrix(-1, ncol=2,nrow=length(result[[run]]$dir2))
	colnames(output[[runname]]$dir2$path) <- c('x','y')
	rownames(output[[runname]]$dir2$path) <- paste('Point',1:nrow(output[[run]]$dir2$path),sep='_')
        for(i in 1:length(result[[run]]$dir2)) {
		output[[runname]]$dir2$path[i,] <- result[[run]]$dir2[[i]]$M
		output[[runname]]$dir2$values[i] <- result[[run]]$dir2[[i]]$val
	}
        
}

output$nrun <- nrun
output$threshold <- Dlim
output$xy <- xy
output$cn <- cn.nb
output$call <- match.call()
class(output) <- 'monmonier'
return(output)
}




##########################
# function plot.monmonier 
##########################
plot.monmonier <- function(x, variable=NULL,displayed.runs=1:x$nrun,
                           add.arrows=TRUE, col='blue',lty=1,bwd=4, clegend=1,csize=0.7,
                           method = c('squaresize','greylevel'),sub='',csub=1,possub='topleft',
                           cneig=1,pixmap=NULL,contour=NULL,area=NULL,add.plot=FALSE,...){

if (!inherits(x, "monmonier")) stop("Use only with 'monmonier' objects")
if(!is.null(variable) & !is.numeric(variable)) stop('If provided, variable must be numeric.\n')
# call <- as.list(x$call) ## no longer used

xy <- x$xy
cpoint <- 0

if(cneig>0) {neig <- nb2neig(x$cn)} else {neig <- NULL}
if(is.null(variable)){
	variable <- rep(1,nrow(xy))
	csize <- 0
	cpoint <- 1
	clegend <- 0
}
s.value(xy,variable,grid=FALSE,include.ori=FALSE,addaxes=FALSE,neig=neig,
        cneig=cneig,clegend=clegend,csize=csize,cpoint=cpoint,pch=20,pixmap=pixmap,
        method=method,sub=sub,csub=csub,possub=possub,add.plot=add.plot)
opar <- par(no.readonly=TRUE)
on.exit(par(mar=opar$mar))
par(mar=c(0,0,0,0))

for(run in displayed.runs){
        obj <- x[[run]]
	if(length(col)!=x$nrun) col <- rep(col,x$nrun)
	if(length(lty)!=x$nrun) lty <- rep(lty,x$nrun)
	if(length(obj$dir1$values) == 0) stop(paste('Monmonier object of run', run, 'is empty (no point in the path)\n'))
	if(length(obj$dir1$values) == 1) points(obj$dir1$path[1],obj$dir1$path[2],pch=20,col=col[run],...)
	else{
	# ceci gere l'épaisseur des lignes
	# la ligne la plus large correspond à la plus grande différence entre deux points
        # comme les deux directions sont gérées séparemment, il en va de meme pour l'epaisseur
        val.1 <- obj$dir1$values
        val.2 <- obj$dir2$values
        n1 <- length(val.1)
        n2 <- length(val.2)
        cex.bwd.1 <- ( val.1[1:(n1-1)] + val.1[2:n1] )/2
        cex.bwd.2 <- ( val.2[1:(n2-1)] + val.2[2:n2] )/2
        cex.bwd.1 <- cex.bwd.1/max(cex.bwd.1)
        cex.bwd.2 <- cex.bwd.2/max(cex.bwd.2)
        
	
		if(add.arrows) {
                  arrows(obj$dir1$path[1:(nrow(obj$dir1$path)-1),1],
                         obj$dir1$path[1:(nrow(obj$dir1$path)-1),2],
                         obj$dir1$path[2:nrow(obj$dir1$path),1],
                         obj$dir1$path[2:nrow(obj$dir1$path),2],
                         lwd=bwd*cex.bwd.1,angle=20,length=0.2,col=col[run],lty=lty[run],...)
                  
                  if(n2>1) arrows(obj$dir2$path[1:(nrow(obj$dir2$path)-1),1],
                                  obj$dir2$path[1:(nrow(obj$dir2$path)-1),2],
                                  obj$dir2$path[2:nrow(obj$dir2$path),1],
                                  obj$dir2$path[2:nrow(obj$dir2$path),2],
                                  lwd=bwd*cex.bwd.2,angle=20,length=0.2,col=col[run],lty=lty[run],...)
                } else {
                  segments(obj$dir1$path[1:(nrow(obj$dir1$path)-1),1],
                           obj$dir1$path[1:(nrow(obj$dir1$path)-1),2],
                           obj$dir1$path[2:nrow(obj$dir1$path),1],
                           obj$dir1$path[2:nrow(obj$dir1$path),2],
                           lwd=bwd*cex.bwd.1,col=col[run],lty=lty[run],...)
                  
                  if(n2>1)segments(obj$dir2$path[1:(nrow(obj$dir2$path)-1),1],
                                   obj$dir2$path[1:(nrow(obj$dir2$path)-1),2],
                                   obj$dir2$path[2:nrow(obj$dir2$path),1],
                                   obj$dir2$path[2:nrow(obj$dir2$path),2],
                                   lwd=bwd*cex.bwd.2,col=col[run],lty=lty[run],...)
                }
	} # end else
} # end for

} # end function



#####################
# print function 
#####################
print.monmonier <- function(x, ...){
cat("\t\n###########################################################")
cat("\t\n# List of paths of maximum differences between neighbours #")
cat("\t\n#           Using a Monmonier based algorithm             #")
cat("\t\n###########################################################\n")
cat('\n$call:')
print(x$call)

cat('\n      # Object content #')
cat("\nClass: ", class(x))
cat('\n$nrun (number of successive runs): ', x$nrun)
if(x$nrun==1)
cat('\n$run1: run of the algorithm')
else if(x$nrun==2)
cat('\n$run1, $run2: runs of the algorithm')
else cat('\n$run1 ... $run',x$nrun, ': runs of the algorithm',sep='')
cat('\n$threshold (minimum difference between neighbours): ', x$threshold)
cat("\n$xy: spatial coordinates")
cat("\n$cn: connection network")

cat('\n\n      # Runs content #')
for(i in 1:x$nrun){
	cat('\n# Run',i)
        # dir 1 #
        cat('\n# First direction')
        cat('\nClass: ', class(x$run1$dir1))
	cat('\n$path:\n')
	print(head(x[[i]]$dir1$path,n=3))
	if(nrow(x[[i]]$dir1$path) >3) cat('...\n')
	cat('\n$values:\n',head(x[[i]]$dir1$values,n=3))
	if(length(x[[i]]$dir1$values)>3) cat(' ...')
        # dir 2 #
        cat('\n# Second direction')
        cat('\nClass: ', class(x$run1$dir2))
	cat('\n$path:\n')
	print(head(x[[i]]$dir2$path,n=3))
	if(nrow(x[[i]]$dir2$path) >3) cat('...\n')
	cat('\n$values:\n',head(x[[i]]$dir2$values,n=3))
	if(length(x[[i]]$dir2$values)>3) cat(' ...')
        
	cat('\n')
	lenTheo <- x$nrun + 5
	if(length(names(x))> lenTheo) {
		cat('Other elements: \n')
		cat(names(x)[(lenTheo+1) : length(x)])
		}
	cat('\n')
	}
}



##############################
# function optimize.monmonier 
##############################
optimize.monmonier <- function(xy,dist,cn,ntry=10,return.best=TRUE, 
                               display.graph=TRUE,threshold=NULL,scanthres=is.null(threshold)){

# a deplacer dans les arguments si on veut permettre l'optimisation sur un objet existant
X <- NULL
#if( any(is.null(xy), is.null(dist),  is.null(cn)) & is.null(X) ) stop("Please provide either xy, dist and cn or a monmonier object (X)")

# si X est un objet monmonier
if(inherits(X,what="monmonier")){
  obj <- as.list(X$call)
  xy <- obj$xy
  dist <- obj$dist
  cn <- obj$cn
}

cn.nb <- cn
cn <- nb2neig(cn)
M <- neig2mat(cn)
D <- as.matrix(dist)
D <- M*D
if(is.null(threshold) || (threshold<0)) {Dlim <- summary(unique(D[D>0]))[5]} else {Dlim <- threshold}

if(scanthres){
  barplot(sort(unique(D)[unique(D) > 0],dec=TRUE),main="Local distances barplot")
  abline(h=Dlim,lty=2)
  mtext("Dashed line indicates present threshold")
  cat("Indicate the threshold (\'d\' for default): ")
  temp <- as.character(readLines(n = 1))
  if(toupper(temp)!="D") { Dlim <- as.numeric(temp) }
}

# série de tests
cat(paste("Boundaries computed (required: ",ntry,")\n",sep=""))

# boucle for obligée pour utiliser aussi pour le cat
tests <- list()
for(i in 0:(ntry-1)){
  tests[[i+1]] <- monmonier(xy, dist, cn.nb,skip=i,scanthres=FALSE,threshold=Dlim)
  cat(paste(1+i," "))
}

bdr.values <- sapply(1:ntry,function(i) sum(c(tests[[i]]$run1$dir1$values,tests[[i]]$run1$dir2$values)) )

# représentation graphique
if(display.graph) barplot(bdr.values,xlab="Local differences skipped",ylab="Sum of all local differences",names.arg=0:(ntry-1))

# retour de la valeur optimale ou de l'objet correspondant
val <- which.max(bdr.values)-1
if(!return.best) {
  return(val)
  } else {
    cat(paste("\nOptimal number of skipped local differences: ",val,"\n"))
    call <- as.list(match.call())
    exp <- bquote( monmonier(xy=.(call$xy),dist=.(call$dist),cn=.(call$cn),skip=.(val),scan=FALSE,thres=.(Dlim)) )
    return(eval(exp))
    }
}
