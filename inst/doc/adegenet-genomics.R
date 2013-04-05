### R code from vignette source 'adegenet-genomics.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: adegenet-genomics.Rnw:105-107
###################################################
library(adegenet)
getClassDef("SNPbin")


###################################################
### code chunk number 2: adegenet-genomics.Rnw:121-122 (eval = FALSE)
###################################################
## new("SNPbin")


###################################################
### code chunk number 3: adegenet-genomics.Rnw:124-125
###################################################
new("SNPbin", multicore=FALSE)


###################################################
### code chunk number 4: adegenet-genomics.Rnw:131-133 (eval = FALSE)
###################################################
## x <- new("SNPbin", c(0,1,1,2,0,0,1))
## x


###################################################
### code chunk number 5: adegenet-genomics.Rnw:135-137
###################################################
x <- new("SNPbin", c(0,1,1,2,0,0,1), multicore=FALSE)
x


###################################################
### code chunk number 6: adegenet-genomics.Rnw:143-146
###################################################
x
ploidy(x) <- 3
x


###################################################
### code chunk number 7: adegenet-genomics.Rnw:150-151
###################################################
x@snp


###################################################
### code chunk number 8: adegenet-genomics.Rnw:154-155
###################################################
as.integer(x)


###################################################
### code chunk number 9: adegenet-genomics.Rnw:162-164
###################################################
dat <- sample(0:1, 1e6, replace=TRUE)
print(object.size(dat),unit="auto")


###################################################
### code chunk number 10: adegenet-genomics.Rnw:166-167 (eval = FALSE)
###################################################
## x <- new("SNPbin", dat, multicore=FALSE)


###################################################
### code chunk number 11: adegenet-genomics.Rnw:169-170
###################################################
load("Robjects/x.snpbin.RData")


###################################################
### code chunk number 12: adegenet-genomics.Rnw:172-173
###################################################
print(object.size(x),unit="auto")


###################################################
### code chunk number 13: adegenet-genomics.Rnw:178-179
###################################################
identical(as.integer(x),dat)


###################################################
### code chunk number 14: adegenet-genomics.Rnw:204-205
###################################################
getClassDef("genlight")


###################################################
### code chunk number 15: adegenet-genomics.Rnw:227-228 (eval = FALSE)
###################################################
## new("genlight")


###################################################
### code chunk number 16: adegenet-genomics.Rnw:230-231
###################################################
new("genlight", multicore=FALSE)


###################################################
### code chunk number 17: adegenet-genomics.Rnw:244-245 (eval = FALSE)
###################################################
## x <- new("genlight", list(indiv1=c(1,1,0,1,1,0), indiv2=c(2,1,1,0,0,0), toto=c(2,2,0,0,4,4)))


###################################################
### code chunk number 18: adegenet-genomics.Rnw:247-248
###################################################
x <- new("genlight", list(indiv1=c(1,1,0,1,1,0), indiv2=c(2,1,1,0,0,0), toto=c(2,2,0,0,4,4)), multicore=FALSE)


###################################################
### code chunk number 19: adegenet-genomics.Rnw:250-252
###################################################
x
ploidy(x)


###################################################
### code chunk number 20: adegenet-genomics.Rnw:257-259
###################################################
as.list(x)
as.matrix(x)


###################################################
### code chunk number 21: adegenet-genomics.Rnw:266-269
###################################################
dat <- lapply(1:50, function(i) sample(c(0,1,NA), 1e5, prob=c(.5, .499, .001), replace=TRUE))
names(dat) <- paste("indiv", 1:length(dat))
print(object.size(dat),unit="auto")


###################################################
### code chunk number 22: adegenet-genomics.Rnw:271-272 (eval = FALSE)
###################################################
## x <- new("genlight", dat)


###################################################
### code chunk number 23: adegenet-genomics.Rnw:274-275
###################################################
load("Robjects/x.genlight.RData")


###################################################
### code chunk number 24: adegenet-genomics.Rnw:277-279
###################################################
print(object.size(x),unit="auto")
object.size(dat)/object.size(x)


###################################################
### code chunk number 25: adegenet-genomics.Rnw:304-305
###################################################
rm(dat)


###################################################
### code chunk number 26: adegenet-genomics.Rnw:343-353
###################################################
dat <- lapply(1:3, function(i) sample(0:2, 10, replace=TRUE))
dat
x <- new("genlight", dat, multicore=FALSE)
x
indNames(x)
indNames(x) <- paste("individual", 1:3)
indNames(x)
locNames(x)
locNames(x) <- paste("SNP",1:nLoc(x),sep=".")
as.matrix(x)


###################################################
### code chunk number 27: adegenet-genomics.Rnw:368-369
###################################################
x


###################################################
### code chunk number 28: adegenet-genomics.Rnw:372-373 (eval = FALSE)
###################################################
## chr(x) <- rep("chr-1", 7)


###################################################
### code chunk number 29: adegenet-genomics.Rnw:376-379
###################################################
chr(x) <- rep("chr-1", 10)
x
chr(x)


###################################################
### code chunk number 30: asmatrix
###################################################
x
as.matrix(x)
as.matrix(x[c(1,3),])
as.matrix(x[, c(TRUE,FALSE)])
as.matrix(x[1:2, c(1,1,1,2,2,2,3,3,3)])


###################################################
### code chunk number 31: seploc
###################################################
x
as.matrix(x)
seploc(x, n.block=2, multicore=FALSE)
lapply(seploc(x, n.block=2, multicore=FALSE),as.matrix)


###################################################
### code chunk number 32: adegenet-genomics.Rnw:417-418
###################################################
lapply(seploc(x, n.block=2, random=TRUE, multicore=FALSE),as.matrix)


###################################################
### code chunk number 33: adegenet-genomics.Rnw:443-444 (eval = FALSE)
###################################################
## file.show(system.file("files/exampleSnpDat.snp",package="adegenet"))


###################################################
### code chunk number 34: readsnp
###################################################
obj <- read.snp(system.file("files/exampleSnpDat.snp",package="adegenet"), chunk=2, multicore=FALSE)
obj
as.matrix(obj, multicore=FALSE)
alleles(obj)
pop(obj)
indNames(obj)


###################################################
### code chunk number 35: adegenet-genomics.Rnw:531-532 (eval = FALSE)
###################################################
## obj <- read.snp("path-to-my-file.snp")


###################################################
### code chunk number 36: adegenet-genomics.Rnw:583-586
###################################################
myPath <- system.file("files/usflu.fasta",package="adegenet")
flu <- fasta2genlight(myPath, chunk=10, multicore=FALSE)
flu


###################################################
### code chunk number 37: adegenet-genomics.Rnw:594-597
###################################################
head(position(flu), 20)
head(alleles(flu), 20)
head(locNames(flu), 20)


###################################################
### code chunk number 38: adegenet-genomics.Rnw:602-607
###################################################
temp <- density(position(flu), bw=10)
plot(temp, type="n", xlab="Position in the alignment", main="Location of the SNPs",
     xlim=c(0,1701))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(position(flu), rep(0, nLoc(flu)), pch="|", col="blue")


###################################################
### code chunk number 39: adegenet-genomics.Rnw:619-621
###################################################
flu <- fasta2genlight(myPath, chunk=10,saveNbAlleles=TRUE, quiet=TRUE, multicore=FALSE)
flu


###################################################
### code chunk number 40: adegenet-genomics.Rnw:626-628
###################################################
head(other(flu)$nb.all.per.loc, 20)
100*mean(unlist(other(flu))>1)


###################################################
### code chunk number 41: adegenet-genomics.Rnw:634-637
###################################################
temp <- table(unlist(other(flu)))
barplot(temp, main="Distribution of the number \nof alleles per loci",
        xlab="Number of alleles", ylab="Number of sites", col=heat.colors(4))


###################################################
### code chunk number 42: adegenet-genomics.Rnw:642-645
###################################################
temp <- temp[-1]
temp <- 100*temp/sum(temp)
round(temp,1)


###################################################
### code chunk number 43: adegenet-genomics.Rnw:685-686
###################################################
glPlot(flu, posi="topleft")


###################################################
### code chunk number 44: adegenet-genomics.Rnw:694-695 (eval = FALSE)
###################################################
## x <- glSim(100, 0, 100, ploidy=2)


###################################################
### code chunk number 45: adegenet-genomics.Rnw:697-698 (eval = FALSE)
###################################################
## plot(x)


###################################################
### code chunk number 46: adegenet-genomics.Rnw:705-706 (eval = FALSE)
###################################################
## x <- glSim(100, 100, ploidy=2, LD=TRUE, block.size=10)


###################################################
### code chunk number 47: adegenet-genomics.Rnw:708-709 (eval = FALSE)
###################################################
## plot(x)


###################################################
### code chunk number 48: adegenet-genomics.Rnw:748-753
###################################################
myFreq <- glMean(flu)
hist(myFreq, proba=TRUE, col="gold", xlab="Allele frequencies",
     main="Distribution of (second) allele frequencies")
temp <- density(myFreq)
lines(temp$x, temp$y*1.8,lwd=3)


###################################################
### code chunk number 49: adegenet-genomics.Rnw:759-765
###################################################
myFreq <- glMean(flu)
myFreq <- c(myFreq, 1-myFreq)
hist(myFreq, proba=TRUE, col="darkseagreen3", xlab="Allele frequencies",
     main="Distribution of allele frequencies", nclass=20)
temp <- density(myFreq, bw=.05)
lines(temp$x, temp$y*2,lwd=3)


###################################################
### code chunk number 50: adegenet-genomics.Rnw:785-791
###################################################
head(glNA(flu),20)
temp <- density(glNA(flu), bw=10)
plot(temp, type="n", xlab="Position in the alignment", main="Location of the missing values (NAs)",
     xlim=c(0,1701))
polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3))
points(glNA(flu), rep(0, nLoc(flu)), pch="|", col="blue")


###################################################
### code chunk number 51: adegenet-genomics.Rnw:819-821
###################################################
x <- glSim(40, 1e4, LD=FALSE, multicore=FALSE)
x


###################################################
### code chunk number 52: adegenet-genomics.Rnw:824-828
###################################################
x <- seploc(x, n.block=10, multicore=FALSE)
class(x)
names(x)
x[1:2]


###################################################
### code chunk number 53: adegenet-genomics.Rnw:832-836
###################################################
lD <- lapply(x, function(e) dist(as.matrix(e)))
class(lD)
names(lD)
class(lD[[1]])


###################################################
### code chunk number 54: adegenet-genomics.Rnw:840-841
###################################################
D <- Reduce("+", lD)


###################################################
### code chunk number 55: adegenet-genomics.Rnw:844-847
###################################################
library(ape)
plot(nj(D), type="fan")
title("A simple NJ tree of simulated genlight data")


###################################################
### code chunk number 56: adegenet-genomics.Rnw:869-873
###################################################
x <- new("genlight", list(a=c(0,0,1,1), b=c(1,1,0,0), c=c(1,1,1,1)), multicore=FALSE)
locNames(x) <- 1:4
x
as.matrix(x)


###################################################
### code chunk number 57: adegenet-genomics.Rnw:877-878
###################################################
glMean(x)


###################################################
### code chunk number 58: adegenet-genomics.Rnw:881-886
###################################################
x <- new("genlight", list(a=c(0,0,2,2), b=c(1,1,0,0), c=c(1,1,1,1)), multicore=FALSE)
locNames(x) <- 1:4
x
as.matrix(x)
ploidy(x)


###################################################
### code chunk number 59: adegenet-genomics.Rnw:903-905
###################################################
M <- as.matrix(x)/ ploidy(x)
apply(M,2,mean)


###################################################
### code chunk number 60: adegenet-genomics.Rnw:913-916
###################################################
as.matrix(x)
glMean(x, alleleAsUnit=TRUE)
glMean(x, alleleAsUnit=FALSE)


###################################################
### code chunk number 61: adegenet-genomics.Rnw:932-933 (eval = FALSE)
###################################################
## pca1 <- glPca(flu)


###################################################
### code chunk number 62: adegenet-genomics.Rnw:935-936
###################################################
load("Robjects/pca1.RData")


###################################################
### code chunk number 63: adegenet-genomics.Rnw:945-946
###################################################
pca1


###################################################
### code chunk number 64: adegenet-genomics.Rnw:953-955
###################################################
scatter(pca1, posi="bottomright")
title("PCA of the US influenza data\n axes 1-2")


###################################################
### code chunk number 65: adegenet-genomics.Rnw:961-966
###################################################
library(ape)
tre <- nj(dist(as.matrix(flu)))
tre
plot(tre, typ="fan", cex=0.7)
title("NJ tree of the US influenza data")


###################################################
### code chunk number 66: adegenet-genomics.Rnw:971-974
###################################################
myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4)
abline(h=0,v=0, col="grey")
add.scatter.eig(pca1$eig[1:40],2,1,2, posi="topright", inset=.05, ratio=.3)


###################################################
### code chunk number 67: adegenet-genomics.Rnw:977-980
###################################################
plot(tre, typ="fan", show.tip=FALSE)
tiplabels(pch=20, col=myCol, cex=4)
title("NJ tree of the US influenza data")


###################################################
### code chunk number 68: adegenet-genomics.Rnw:1003-1005 (eval = FALSE)
###################################################
## x <- glSim(100, 1e4, 50)
## dapc1 <- dapc(x, n.pca=10, n.da=1)


###################################################
### code chunk number 69: adegenet-genomics.Rnw:1008-1009
###################################################
load("Robjects/x.dapc1.RData")


###################################################
### code chunk number 70: adegenet-genomics.Rnw:1015-1017
###################################################
scatter(dapc1,scree.da=FALSE, bg="white", posi.pca="topright", legend=TRUE,
        txt.leg=paste("group", 1:2), col=c("red","blue"))


###################################################
### code chunk number 71: adegenet-genomics.Rnw:1021-1022
###################################################
compoplot(dapc1, col=c("red","blue"),lab="", txt.leg=paste("group", 1:2), ncol=2)


###################################################
### code chunk number 72: adegenet-genomics.Rnw:1026-1027 (eval = FALSE)
###################################################
## loadingplot(dapc1$var.contr, thres=1e-3)


###################################################
### code chunk number 73: adegenet-genomics.Rnw:1036-1037 (eval = FALSE)
###################################################
## loadingplot(tail(dapc1$var.contr[,1],100), thres=1e-3)


