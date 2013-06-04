### R code from vignette source 'adegenet-spca.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: adegenet-spca.Rnw:157-160
###################################################
library(ade4)
library(adehabitat)
library(spdep)


###################################################
### code chunk number 2: adegenet-spca.Rnw:162-169
###################################################
library(ade4)
library(adegenet)
library(adehabitat)
data(spcaIllus)
obj <- spcaIllus$dat2A
obj
head(truenames(obj[loc="L01"])$tab)


###################################################
### code chunk number 3: adegenet-spca.Rnw:179-180
###################################################
args(spca)


###################################################
### code chunk number 4: adegenet-spca.Rnw:197-198
###################################################
mySpca <- spca(obj, ask=FALSE, type=1, scannf=FALSE)


###################################################
### code chunk number 5: adegenet-spca.Rnw:218-219 (eval = FALSE)
###################################################
## mySpca <- spca(obj,type=1,ask=FALSE,scannf=FALSE)


###################################################
### code chunk number 6: adegenet-spca.Rnw:224-225
###################################################
mySpca <- spca(obj,type=5,d1=0,d2=2,scannf=FALSE)


###################################################
### code chunk number 7: adegenet-spca.Rnw:237-241
###################################################
myCn <- chooseCN(obj$other$xy, type=6, k=10, plot=FALSE)
myCn
class(myCn)
mySpca2 <- spca(obj,cn=myCn,scannf=FALSE)


###################################################
### code chunk number 8: adegenet-spca.Rnw:253-254
###################################################
barplot(mySpca$eig,main="Eigenvalues of sPCA", col=rep(c("red","grey"),c(1,100)))


###################################################
### code chunk number 9: adegenet-spca.Rnw:278-281
###################################################
mySpca <- spca(obj,type=1,scannf=FALSE,plot.nb=FALSE,nfposi=1,nfnega=0)
class(mySpca)
mySpca


###################################################
### code chunk number 10: adegenet-spca.Rnw:293-300
###################################################
head(mySpca$eig)
tail(mySpca$eig)
length(mySpca$eig)
myPal <- colorRampPalette(c("red","grey","blue"))
barplot(mySpca$eig, main="A variant of the plot\n of sPCA eigenvalues", col=myPal(length(mySpca$eig)))
legend("topright", fill=c("red","blue"), leg=c("Global structures", "Local structures"))
abline(h=0,col="grey")


###################################################
### code chunk number 11: adegenet-spca.Rnw:306-309
###################################################
head(mySpca$c1)
tail(mySpca$c1)
dim(mySpca$c1)


###################################################
### code chunk number 12: adegenet-spca.Rnw:315-318
###################################################
head(mySpca$li)
tail(mySpca$li)
dim(mySpca$li)


###################################################
### code chunk number 13: adegenet-spca.Rnw:324-327
###################################################
head(mySpca$ls)
tail(mySpca$ls)
dim(mySpca$ls)


###################################################
### code chunk number 14: adegenet-spca.Rnw:337-338
###################################################
mySpca$as


###################################################
### code chunk number 15: screeplot
###################################################
screeplot(mySpca)


###################################################
### code chunk number 16: globalrtest
###################################################
myGtest <- global.rtest(obj$tab,mySpca$lw,nperm=99)
myGtest
plot(myGtest)


###################################################
### code chunk number 17: localrtest
###################################################
myLtest <- local.rtest(obj$tab,mySpca$lw,nperm=99)
myLtest
plot(myLtest)


###################################################
### code chunk number 18: plotspca
###################################################
plot(mySpca)


###################################################
### code chunk number 19: colorplot
###################################################
colorplot(mySpca,cex=3,main="colorplot of mySpca, first global score")


###################################################
### code chunk number 20: adegenet-spca.Rnw:480-483
###################################################
library(akima)
x <- other(obj)$xy[,1]
y <- other(obj)$xy[,2]


###################################################
### code chunk number 21: adegenet-spca.Rnw:485-487
###################################################
temp <- interp(x, y, mySpca$li[,1])
image(temp)


###################################################
### code chunk number 22: adegenet-spca.Rnw:493-497
###################################################
interpX <- seq(min(x),max(x),le=200)
interpY <- seq(min(y),max(y),le=200)
temp <- interp(x, y, mySpca$ls[,1], xo=interpX, yo=interpY)
image(temp)


###################################################
### code chunk number 23: adegenet-spca.Rnw:502-508
###################################################
myPal <- colorRampPalette(c("firebrick2", "white", "lightslateblue"))
annot <- function(){
    title("sPCA - interpolated map of individual scores")
    points(x,y)
}
filled.contour(temp, color.pal=myPal, nlev=50, key.title=title("lagged \nscore 1"), plot.title=annot())


###################################################
### code chunk number 24: adegenet-spca.Rnw:519-524
###################################################
myLoadings <- mySpca$c1[,1]^2
names(myLoadings) <- rownames(mySpca$c1)
loadingplot(myLoadings, xlab="Alleles",
            ylab="Weight of the alleles",
            main="Contribution of alleles \n to the first sPCA axis")


###################################################
### code chunk number 25: adegenet-spca.Rnw:538-543
###################################################
temp <- loadingplot(myLoadings, threshold=quantile(myLoadings, 0.95),
                    xlab="Alleles",ylab="Weight of the alleles",
                    main="Contribution of alleles \n to the first sPCA axis",
                    fac=obj$loc.fac, cex.fac=0.6)
temp


###################################################
### code chunk number 26: boxplot
###################################################
boxplot(myLoadings~obj$loc.fac, las=3, ylab="Contribution", xlab="Marker",
        main="Contributions by markers \nto the first global score", col="grey")


###################################################
### code chunk number 27: adegenet-spca.Rnw:589-591
###################################################
data(rupica)
rupica


###################################################
### code chunk number 28: adegenet-spca.Rnw:600-602
###################################################
rupica$other$showBauges()
points(rupica$other$xy, col="red",pch=20)


###################################################
### code chunk number 29: adegenet-spca.Rnw:610-613
###################################################
rupica$other$showBauges()
s.kde2d(rupica$other$xy,add.plot=TRUE)
points(rupica$other$xy, col="red",pch=20)


###################################################
### code chunk number 30: adegenet-spca.Rnw:632-635
###################################################
rupica.smry <- summary(rupica)
plot(rupica.smry$Hobs, rupica.smry$Hexp, main="Observed vs expected heterozygosity")
abline(0,1,col="red")


###################################################
### code chunk number 31: adegenet-spca.Rnw:641-642
###################################################
t.test(rupica.smry$Hexp, rupica.smry$Hobs,paired=TRUE,var.equal=TRUE)


###################################################
### code chunk number 32: adegenet-spca.Rnw:655-659
###################################################
rupica.X <- scaleGen(rupica)
rupica.pca1 <- dudi.pca(rupica.X, cent=FALSE, scale=FALSE, scannf=FALSE, nf=2)
barplot(rupica.pca1$eig, main="Rupica dataset - PCA eigenvalues",
        col=heat.colors(length(rupica.pca1$eig)))


###################################################
### code chunk number 33: adegenet-spca.Rnw:667-668
###################################################
rupica.pca1


###################################################
### code chunk number 34: adegenet-spca.Rnw:682-685
###################################################
s.label(rupica.pca1$li)
s.kde2d(rupica.pca1$li, add.p=TRUE, cpoint=0)
add.scatter.eig(rupica.pca1$eig,2,1,2)


###################################################
### code chunk number 35: adegenet-spca.Rnw:690-691
###################################################
loadingplot(rupica.pca1$c1^2)


###################################################
### code chunk number 36: adegenet-spca.Rnw:698-703
###################################################
X <- truenames(rupica)
class(X)
dim(X)
bm203.221 <- X[,"Bm203.221"]
table(bm203.221)


###################################################
### code chunk number 37: adegenet-spca.Rnw:708-709
###################################################
rownames(X)[bm203.221==0.5]


###################################################
### code chunk number 38: svaluedem
###################################################
s.value(cbind(1:11,rep(1,11)), -5:5, cleg=0)
text(1:11,rep(1,11), -5:5, col="red",cex=1.5)


###################################################
### code chunk number 39: adegenet-spca.Rnw:735-739
###################################################
showBauges <- rupica$other$showBauges
showBauges()
s.value(rupica$other$xy, rupica.pca1$li[,1], add.p=TRUE, cleg=0.5)
title("PCA - first PC",col.main="yellow" ,line=-2, cex.main=2)


###################################################
### code chunk number 40: adegenet-spca.Rnw:742-745
###################################################
showBauges()
s.value(rupica$other$xy, rupica.pca1$li[,2], add.p=TRUE, csize=0.7)
title("PCA - second PC",col.main="yellow" ,line=-2, cex.main=2)


###################################################
### code chunk number 41: adegenet-spca.Rnw:757-758
###################################################
rupica.graph <- chooseCN(rupica$other$xy,type=5,d1=0,d2=2300, plot=FALSE, res="listw")


###################################################
### code chunk number 42: adegenet-spca.Rnw:761-764
###################################################
rupica.graph
plot(rupica.graph, rupica$other$xy)
title("rupica.graph")


###################################################
### code chunk number 43: adegenet-spca.Rnw:768-770
###################################################
pc1.mctest <- moran.mc(rupica.pca1$li[,1], rupica.graph, 999)
plot(pc1.mctest)


###################################################
### code chunk number 44: adegenet-spca.Rnw:776-777
###################################################
moran.plot(rupica.pca1$li[,1], rupica.graph)


###################################################
### code chunk number 45: adegenet-spca.Rnw:786-788
###################################################
pc2.mctest <- moran.mc(rupica.pca1$li[,2], rupica.graph, 999)
plot(pc2.mctest)


###################################################
### code chunk number 46: adegenet-spca.Rnw:808-810
###################################################
mtest <- mantel.randtest(dist(rupica.X), dist(rupica$other$xy))
plot(mtest, nclass=30)


###################################################
### code chunk number 47: adegenet-spca.Rnw:826-828
###################################################
rupica.spca1 <- spca(rupica, cn=rupica.graph,scannf=FALSE, nfposi=2,nfnega=0)
barplot(rupica.spca1$eig, col=rep(c("red","grey"), c(2,1000)), main="rupica dataset - sPCA eigenvalues")


###################################################
### code chunk number 48: adegenet-spca.Rnw:835-836
###################################################
rupica.spca1


###################################################
### code chunk number 49: adegenet-spca.Rnw:845-846
###################################################
screeplot(rupica.spca1)


###################################################
### code chunk number 50: adegenet-spca.Rnw:859-862
###################################################
showBauges()
s.value(rupica$other$xy, rupica.spca1$ls[,1], add.p=TRUE, csize=0.7)
title("sPCA - first PC",col.main="yellow" ,line=-2, cex.main=2)


###################################################
### code chunk number 51: adegenet-spca.Rnw:870-873
###################################################
showBauges()
s.value(rupica$other$xy, rupica.spca1$ls[,2], add.p=TRUE, csize=0.7)
title("sPCA - second PC",col.main="yellow" ,line=-2, cex.main=2)


###################################################
### code chunk number 52: adegenet-spca.Rnw:884-887
###################################################
showBauges()
colorplot(rupica$other$xy, rupica.spca1$ls, axes=1:2, transp=TRUE, add=TRUE, cex=3)
title("sPCA - colorplot of PC 1 and 2\n(lagged scores)", col.main="yellow", line=-2, cex=2)


