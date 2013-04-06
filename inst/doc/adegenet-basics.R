### R code from vignette source 'adegenet-basics.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: adegenet-basics.Rnw:120-121
###################################################
R.version.string


###################################################
### code chunk number 2: adegenet-basics.Rnw:125-126 (eval = FALSE)
###################################################
## install.packages("adegenet", dep=TRUE)


###################################################
### code chunk number 3: adegenet-basics.Rnw:130-133
###################################################
library(ape)
library(seqinr)
library(genetics)


###################################################
### code chunk number 4: adegenet-basics.Rnw:135-136
###################################################
library(adegenet)


###################################################
### code chunk number 5: adegenet-basics.Rnw:140-141
###################################################
packageDescription("adegenet", fields = "Version")


###################################################
### code chunk number 6: adegenet-basics.Rnw:155-156 (eval = FALSE)
###################################################
## help.search("Hardy-Weinberg")


###################################################
### code chunk number 7: adegenet-basics.Rnw:163-164 (eval = FALSE)
###################################################
## ?spca


###################################################
### code chunk number 8: adegenet-basics.Rnw:182-183 (eval = FALSE)
###################################################
## adegenetWeb()


###################################################
### code chunk number 9: adegenet-basics.Rnw:190-191 (eval = FALSE)
###################################################
## ?adegenet


###################################################
### code chunk number 10: adegenet-basics.Rnw:195-196 (eval = FALSE)
###################################################
## help.start()


###################################################
### code chunk number 11: genind
###################################################
data(nancycats)
is.genind(nancycats)
nancycats


###################################################
### code chunk number 12: adegenet-basics.Rnw:256-257
###################################################
getClassDef("genind")


###################################################
### code chunk number 13: adegenet-basics.Rnw:287-288
###################################################
nancycats$tab[10:18,1:10]


###################################################
### code chunk number 14: adegenet-basics.Rnw:295-296
###################################################
nancycats$loc.names


###################################################
### code chunk number 15: adegenet-basics.Rnw:299-300
###################################################
nancycats$all.names[[3]]


###################################################
### code chunk number 16: adegenet-basics.Rnw:325-328
###################################################
head(pop(nancycats))
catpop <- genind2genpop(nancycats)
catpop


###################################################
### code chunk number 17: adegenet-basics.Rnw:341-345
###################################################
obj <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"))
obj$call
toto <- eval(obj$call)
identical(obj,toto)


###################################################
### code chunk number 18: adegenet-basics.Rnw:353-356
###################################################
catpop
is.genpop(catpop)
catpop$tab[1:5,1:10]


###################################################
### code chunk number 19: adegenet-basics.Rnw:362-363
###################################################
getClassDef("genpop")


###################################################
### code chunk number 20: adegenet-basics.Rnw:393-396
###################################################
head(indNames(nancycats),10)
indNames(nancycats) <- paste("cat", 1:nInd(nancycats),sep=".")
head(indNames(nancycats),10)


###################################################
### code chunk number 21: adegenet-basics.Rnw:400-401
###################################################
locNames(nancycats)


###################################################
### code chunk number 22: adegenet-basics.Rnw:404-406
###################################################
temp <- locNames(nancycats, withAlleles=TRUE)
head(temp, 10)


###################################################
### code chunk number 23: adegenet-basics.Rnw:412-416
###################################################
obj <- nancycats[sample(1:50,10)]
pop(obj)
pop(obj) <- rep("newPop",10)
pop(obj)


###################################################
### code chunk number 24: adegenet-basics.Rnw:421-422 (eval = FALSE)
###################################################
## obj@pop <- rep("newPop",10)


###################################################
### code chunk number 25: import
###################################################
obj1 <- read.genetix(system.file("files/nancycats.gtx",package="adegenet"))
obj2 <- import2genind(system.file("files/nancycats.gtx", package="adegenet"))
all.equal(obj1,obj2)



###################################################
### code chunk number 26: adegenet-basics.Rnw:490-496
###################################################
temp <- lapply(1:30, function(i) sample(1:9, 4, replace=TRUE))
temp <- sapply(temp, paste, collapse="")
temp <- matrix(temp, nrow=10, dimnames=list(paste("ind",1:10), paste("loc",1:3)))
temp
obj <- df2genind(temp, ploidy=4, sep="")
obj


###################################################
### code chunk number 27: adegenet-basics.Rnw:503-504
###################################################
genind2df(obj, sep="|")


###################################################
### code chunk number 28: aflpread
###################################################
dat <- read.table(system.file("files/AFLP.txt",package="adegenet"), header=TRUE)
dat


###################################################
### code chunk number 29: adegenet-basics.Rnw:546-549
###################################################
obj <- genind(dat, ploidy=1, type="PA")
obj
truenames(obj)


###################################################
### code chunk number 30: adegenet-basics.Rnw:553-555
###################################################
pop(obj) <- rep(c('a','b'),4:3)
summary(obj)


###################################################
### code chunk number 31: adegenet-basics.Rnw:560-563
###################################################
obj2 <- genind2genpop(obj)
obj2
truenames(obj2)


###################################################
### code chunk number 32: adegenet-basics.Rnw:568-570
###################################################
objNoNa <- na.replace(obj,met=0)
objNoNa@tab


###################################################
### code chunk number 33: pcaaflp
###################################################
library(ade4)
pca1 <- dudi.pca(objNoNa,scannf=FALSE,scale=FALSE)
temp <- as.integer(pop(objNoNa))
myCol <- transp(c("blue","red"),.7)[temp]
myPch <- c(15,17)[temp]
plot(pca1$li, col=myCol, cex=3, pch=myPch)
abline(h=0,v=0,col="grey",lty=2)
s.arrow(pca1$c1, add.plot=TRUE)
legend("topright", pch=c(15,17), col=transp(c("blue","red"),.7), leg=c("Group A","Group B"), pt.cex=2)


###################################################
### code chunk number 34: adegenet-basics.Rnw:622-628
###################################################
dat <- matrix(sample(c("a","t","g","c"), 15, replace=TRUE),nrow=3)
rownames(dat) <- paste("genot.", 1:3)
colnames(dat) <- 1:5
dat
obj <- df2genind(dat, ploidy=1)
truenames(obj)


###################################################
### code chunk number 35: adegenet-basics.Rnw:655-661
###################################################
library(ape)
ref <- c("U15717", "U15718", "U15719", "U15720",
         "U15721", "U15722", "U15723", "U15724")
myDNA <- read.GenBank(ref)
myDNA
class(myDNA)


###################################################
### code chunk number 36: adegenet-basics.Rnw:669-671
###################################################
obj <- DNAbin2genind(myDNA, polyThres=0.01)
obj


###################################################
### code chunk number 37: adegenet-basics.Rnw:676-677
###################################################
head(locNames(obj))


###################################################
### code chunk number 38: seqinr1
###################################################
library(seqinr)
mase.res   <- read.alignment(file = system.file("sequences/test.mase",package = "seqinr"), format = "mase")
mase.res
x <- alignment2genind(mase.res)
x


###################################################
### code chunk number 39: adegenet-basics.Rnw:712-713
###################################################
head(locNames(x))


###################################################
### code chunk number 40: adegenet-basics.Rnw:716-719
###################################################
tabAA <- genind2df(x)
dim(tabAA)
tabAA[, 1:20]


###################################################
### code chunk number 41: adegenet-basics.Rnw:722-723
###################################################
table(unlist(tabAA))


###################################################
### code chunk number 42: adegenet-basics.Rnw:733-735
###################################################
D <- dist(truenames(x))
D


###################################################
### code chunk number 43: njAA
###################################################
library(ape)
tre <- nj(D)
par(xpd=TRUE)
plot(tre, type="unrooted", edge.w=2)
edgelabels(tex=round(tre$edge.length,1), bg=rgb(.8,.8,1,.8))


###################################################
### code chunk number 44: adegenet-basics.Rnw:751-754
###################################################
pco1 <- dudi.pco(D, scannf=FALSE,nf=2)
scatter(pco1, posi="bottomright")
title("Principal Coordinate Analysis\n-based on proteic distances-")


###################################################
### code chunk number 45: adegenet-basics.Rnw:776-779
###################################################
library(ade4)
data(microsatt)
microsatt$tab[10:15,12:15]


###################################################
### code chunk number 46: adegenet-basics.Rnw:784-787
###################################################
toto <- genpop(microsatt$tab)
toto
summary(toto)


###################################################
### code chunk number 47: genind2genotype
###################################################
obj <- genind2genotype(nancycats)
class(obj)
obj[1:4,1:5]
class(obj$fca8)


###################################################
### code chunk number 48: genind2hierfstat
###################################################
obj <- genind2hierfstat(nancycats)
class(obj)
obj[1:4,1:5]


###################################################
### code chunk number 49: genind2df
###################################################
obj <- genind2df(nancycats)
obj[1:5,1:5]


###################################################
### code chunk number 50: adegenet-basics.Rnw:841-842
###################################################
genind2df(nancycats,sep="|")[1:5,1:5]


###################################################
### code chunk number 51: adegenet-basics.Rnw:867-873
###################################################
data(microbov)
toto <- genind2genpop(microbov)
toto
toto@pop.names
titi <- toto[1:3,]
titi@pop.names


###################################################
### code chunk number 52: adegenet-basics.Rnw:880-882
###################################################
tata <- titi[,loc="L03"]
tata


###################################################
### code chunk number 53: seploc
###################################################
data(nancycats)
sepCats <- seploc(nancycats)
class(sepCats)
names(sepCats)
sepCats$fca45


###################################################
### code chunk number 54: seppop
###################################################
data(microbov)
obj <- seppop(microbov)
class(obj)
names(obj)
obj$Borgou


###################################################
### code chunk number 55: sepultim
###################################################
obj <- lapply(obj,seploc)
names(obj)
class(obj$Borgou)
names(obj$Borgou)
obj$Borgou$INRA63


###################################################
### code chunk number 56: repool
###################################################
obj <- seppop(microbov)
names(obj)
newObj <- repool(obj$Borgou, obj$Charolais)
newObj
newObj$pop.names


###################################################
### code chunk number 57: sumry
###################################################
toto <- summary(nancycats)
names(toto)

par(mfrow=c(2,2))

plot(toto$pop.eff,toto$pop.nall,xlab="Colonies sample size",ylab="Number of alleles",main="Alleles numbers and sample sizes",type="n")
text(toto$pop.eff,toto$pop.nall,lab=names(toto$pop.eff))

barplot(toto$loc.nall,ylab="Number of alleles", main="Number of alleles per locus")

barplot(toto$Hexp-toto$Hobs,main="Heterozygosity: expected-observed",ylab="Hexp - Hobs")

barplot(toto$pop.eff,main="Sample sizes per population",ylab="Number of genotypes",las=3)


###################################################
### code chunk number 58: adegenet-basics.Rnw:983-985
###################################################
bartlett.test(list(toto$Hexp,toto$Hobs))
t.test(toto$Hexp,toto$Hobs,pair=T,var.equal=TRUE,alter="greater")


###################################################
### code chunk number 59: HWE
###################################################
toto <- HWE.test.genind(nancycats,res="matrix")
dim(toto)


###################################################
### code chunk number 60: adegenet-basics.Rnw:1006-1009
###################################################
colnames(toto)
idx <- which(toto<0.0001,TRUE)
idx


###################################################
### code chunk number 61: adegenet-basics.Rnw:1014-1016
###################################################
toto <- HWE.test.genind(nancycats,res="full")
mapply(function(i,j) toto[[i]][[j]], idx[,2], idx[,1], SIMPLIFY=FALSE)


###################################################
### code chunk number 62: adegenet-basics.Rnw:1087-1090
###################################################
data(nancycats)
matFst <- pairwise.fst(nancycats[1:50, treatOther=FALSE])
matFst


###################################################
### code chunk number 63: adegenet-basics.Rnw:1094-1095
###################################################
is.euclid(matFst)


###################################################
### code chunk number 64: adegenet-basics.Rnw:1144-1147
###################################################
data(microbov)
sal <- seppop(microbov)$Salers
sal


###################################################
### code chunk number 65: adegenet-basics.Rnw:1150-1154
###################################################
temp <- inbreeding(sal, N=100)
class(temp)
head(names(temp))
head(temp[[1]],20)


###################################################
### code chunk number 66: adegenet-basics.Rnw:1158-1159
###################################################
Fbar <- sapply(temp, mean)


###################################################
### code chunk number 67: adegenet-basics.Rnw:1161-1162
###################################################
hist(Fbar, col="firebrick", main="Average inbreeding in Salers cattles")


###################################################
### code chunk number 68: adegenet-basics.Rnw:1167-1170
###################################################
which(Fbar>0.4)
F <- inbreeding(sal, res.type="function")[which(Fbar>0.4)]
F


###################################################
### code chunk number 69: adegenet-basics.Rnw:1174-1175
###################################################
plot(F$FRBTSAL9266, main=paste("Inbreeding of individual",names(F)), xlab="Inbreeding (F)", ylab="Probability density")


###################################################
### code chunk number 70: pcaexpl
###################################################
data(microbov)
sum(is.na(microbov$tab))


###################################################
### code chunk number 71: adegenet-basics.Rnw:1296-1300
###################################################
X <- scaleGen(microbov, missing="mean")
class(X)
dim(X)
X[1:5,1:5]


###################################################
### code chunk number 72: adegenet-basics.Rnw:1310-1312
###################################################
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=3)
barplot(pca1$eig[1:50],main="PCA eigenvalues", col=heat.colors(50))


###################################################
### code chunk number 73: adegenet-basics.Rnw:1314-1315
###################################################
pca1


###################################################
### code chunk number 74: adegenet-basics.Rnw:1328-1331
###################################################
s.label(pca1$li)
title("PCA of microbov dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)


###################################################
### code chunk number 75: adegenet-basics.Rnw:1336-1339
###################################################
s.class(pca1$li, pop(microbov))
title("PCA of microbov dataset\naxes 1-2")
add.scatter.eig(pca1$eig[1:20], 3,1,2)


###################################################
### code chunk number 76: adegenet-basics.Rnw:1344-1347
###################################################
s.class(pca1$li,pop(microbov),xax=1,yax=3,sub="PCA 1-3",csub=2)
title("PCA of microbov dataset\naxes 1-3")
add.scatter.eig(pca1$eig[1:20],nf=3,xax=1,yax=3)


###################################################
### code chunk number 77: adegenet-basics.Rnw:1356-1358
###################################################
col <- rainbow(length(levels(pop(microbov))))
s.class(pca1$li, pop(microbov),xax=1,yax=3, col=transp(col,.6), axesell=FALSE, cstar=0, cpoint=3, grid=FALSE)


###################################################
### code chunk number 78: adegenet-basics.Rnw:1367-1370
###################################################
colorplot(pca1$li, pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 2")
title("PCA of microbov dataset\naxes 1-2")
abline(v=0,h=0,col="grey", lty=2)


###################################################
### code chunk number 79: adegenet-basics.Rnw:1380-1383
###################################################
colorplot(pca1$li[c(1,3)], pca1$li, transp=TRUE, cex=3, xlab="PC 1", ylab="PC 3")
title("PCA of microbov dataset\naxes 1-3")
abline(v=0,h=0,col="grey", lty=2)


###################################################
### code chunk number 80: caexpl
###################################################
data(microbov)
obj <- genind2genpop(microbov,missing="chi2")
ca1 <- dudi.coa(as.data.frame(obj$tab),scannf=FALSE,nf=3)
barplot(ca1$eig,main="Correspondance Analysis eigenvalues", col=heat.colors(length(ca1$eig)))


###################################################
### code chunk number 81: adegenet-basics.Rnw:1406-1408
###################################################
s.label(ca1$li,lab=obj$pop.names,sub="CA 1-2",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=1,yax=2,posi="bottomright")


###################################################
### code chunk number 82: adegenet-basics.Rnw:1411-1413
###################################################
s.label(ca1$li,xax=1,yax=3,lab=obj$pop.names,sub="CA 1-3",csub=2)
add.scatter.eig(ca1$eig,nf=3,xax=2,yax=3,posi="topleft")


###################################################
### code chunk number 83: ibd
###################################################
data(nancycats)
toto <- genind2genpop(nancycats,miss="0")
Dgen <- dist.genpop(toto,method=2)
Dgeo <- dist(nancycats$other$xy)
ibd <- mantel.randtest(Dgen,Dgeo)
ibd


###################################################
### code chunk number 84: adegenet-basics.Rnw:1504-1505
###################################################
plot(ibd)


###################################################
### code chunk number 85: adegenet-basics.Rnw:1517-1523
###################################################
data(spcaIllus)
x <- spcaIllus$dat2B
Dgen <- dist(x$tab)
Dgeo <- dist(other(x)$xy)
ibd <- mantel.randtest(Dgen,Dgeo)
ibd


###################################################
### code chunk number 86: adegenet-basics.Rnw:1525-1526
###################################################
plot(ibd)


###################################################
### code chunk number 87: adegenet-basics.Rnw:1541-1543
###################################################
plot(Dgeo, Dgen)
abline(lm(Dgen~Dgeo), col="red",lty=2)


###################################################
### code chunk number 88: adegenet-basics.Rnw:1552-1558 (eval = FALSE)
###################################################
## dens <- kde2d(Dgeo,Dgen, n=300, lims=c(-.1, 1.5,-.5,4))
## myPal <- colorRampPalette(c("white","blue","gold", "orange", "red"))
## plot(Dgeo, Dgen, pch=20,cex=.5)
## image(dens, col=transp(myPal(300),.7), add=TRUE)
## abline(lm(Dgen~Dgeo))
## title("Isolation by distance plot")


###################################################
### code chunk number 89: mon1
###################################################
data(sim2pop)
sim2pop
summary(sim2pop$pop)

temp <- sim2pop$pop
levels(temp) <- c(3,5)
temp <- as.numeric(as.character(temp))
plot(sim2pop$other$xy,pch=temp,cex=1.5,xlab='x',ylab='y')
legend("topright",leg=c("Pop A", "Pop B"),pch=c(3,5))


###################################################
### code chunk number 90: mon2
###################################################
args(monmonier)


###################################################
### code chunk number 91: mon3
###################################################
D <- dist(sim2pop$tab)


###################################################
### code chunk number 92: mon4
###################################################
gab <- chooseCN(sim2pop$other$xy,ask=FALSE,type=2)


###################################################
### code chunk number 93: mon5 (eval = FALSE)
###################################################
## mon1 <- monmonier(sim2pop$other$xy,D,gab)


###################################################
### code chunk number 94: adegenet-basics.Rnw:1644-1645
###################################################
pairwise.fst(sim2pop)


###################################################
### code chunk number 95: adegenet-basics.Rnw:1652-1653
###################################################
replicate(10, pairwise.fst(sim2pop, pop=sample(pop(sim2pop))))


###################################################
### code chunk number 96: mon6
###################################################
library(ade4)
pco1 <- dudi.pco(D,scannf=FALSE,nf=1)
barplot(pco1$eig,main="Eigenvalues")


###################################################
### code chunk number 97: mon7
###################################################
D <- dist(pco1$li)


###################################################
### code chunk number 98: mon8 (eval = FALSE)
###################################################
## mon1 <- monmonier(sim2pop$other$xy,D,gab)


###################################################
### code chunk number 99: mon9
###################################################
mon1 <- monmonier(sim2pop$other$xy,D,gab,scan=FALSE)
mon1


###################################################
### code chunk number 100: mon10
###################################################
names(mon1)
names(mon1$run1)
mon1$run1$dir1


###################################################
### code chunk number 101: adegenet-basics.Rnw:1700-1701
###################################################
coords.monmonier(mon1)


###################################################
### code chunk number 102: adegenet-basics.Rnw:1710-1711
###################################################
plot(mon1)


###################################################
### code chunk number 103: adegenet-basics.Rnw:1716-1722
###################################################
plot(mon1,add.arrows=FALSE,bwd=8)
temp <- sim2pop$pop
levels(temp) <- c(3,5)
temp <- as.numeric(as.character(temp))
points(sim2pop$other$xy,pch=temp,cex=1.3)
legend("topright",leg=c("Pop A", "Pop B"),pch=c(3,5))


###################################################
### code chunk number 104: hybr
###################################################
temp <- seppop(microbov)
names(temp)
salers <- temp$Salers
zebu <- temp$Zebu
zebler <- hybridize(salers, zebu, n=40, pop="zebler")


###################################################
### code chunk number 105: adegenet-basics.Rnw:1753-1756
###################################################
F2 <- hybridize(salers, zebler, n=40)
F3 <- hybridize(salers, F2, n=40)
F4 <- hybridize(salers, F3, n=40)


