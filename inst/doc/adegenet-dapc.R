### R code from vignette source 'adegenet-dapc.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: adegenet-dapc.Rnw:151-155
###################################################
library(adegenet)
data(dapcIllus)
class(dapcIllus)
names(dapcIllus)


###################################################
### code chunk number 2: adegenet-dapc.Rnw:159-161
###################################################
x <- dapcIllus$a
x


###################################################
### code chunk number 3: adegenet-dapc.Rnw:167-168
###################################################
grp <- find.clusters(x, n.pca=200, n.clust=6)


###################################################
### code chunk number 4: adegenet-dapc.Rnw:170-171 (eval = FALSE)
###################################################
## grp <- find.clusters(x, max.n.clust=40)


###################################################
### code chunk number 5: adegenet-dapc.Rnw:196-201
###################################################
names(grp)
head(grp$Kstat, 8)
grp$stat
head(grp$grp, 10)
grp$size


###################################################
### code chunk number 6: adegenet-dapc.Rnw:209-211
###################################################
table(pop(x), grp$grp)
table.value(table(pop(x), grp$grp), col.lab=paste("inf", 1:6), row.lab=paste("ori", 1:6))


###################################################
### code chunk number 7: adegenet-dapc.Rnw:296-297
###################################################
dapc1 <- dapc(x, grp$grp, n.pca=40, n.da=100)


###################################################
### code chunk number 8: adegenet-dapc.Rnw:299-300 (eval = FALSE)
###################################################
## dapc1 <- dapc(x, grp$grp)


###################################################
### code chunk number 9: adegenet-dapc.Rnw:329-330
###################################################
dapc1


###################################################
### code chunk number 10: adegenet-dapc.Rnw:340-341
###################################################
scatter(dapc1)


###################################################
### code chunk number 11: adegenet-dapc.Rnw:369-370
###################################################
scatter(dapc1, posi.da="bottomright", bg="white", pch=17:22)


###################################################
### code chunk number 12: adegenet-dapc.Rnw:375-377
###################################################
myCol <- c("darkblue","purple","green","orange","red","blue")
scatter(dapc1, posi.da="bottomright",  bg="white", pch=17:22, cstar=0, col=myCol, scree.pca=TRUE, posi.pca="bottomleft")


###################################################
### code chunk number 13: adegenet-dapc.Rnw:383-385
###################################################
scatter(dapc1, scree.da=FALSE, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4,
        cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:6))


###################################################
### code chunk number 14: adegenet-dapc.Rnw:396-412
###################################################
scatter(dapc1, ratio.pca=0.3, bg="white", pch=20,  cell=0, cstar=0, col=myCol, solid=.4,
        cex=3, clab=0, mstree=TRUE, scree.da=FALSE,
        posi.pca="bottomright", leg=TRUE, txt.leg=paste("Cluster",1:6))

par(xpd=TRUE)
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=8, col="black")
points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4, cex=3, lwd=2, col=myCol)

myInset <- function(){
    temp <- dapc1$pca.eig
    temp <- 100* cumsum(temp)/sum(temp)
    plot(temp, col=rep(c("black","lightgrey"), c(dapc1$n.pca,1000)), ylim=c(0,100),
         xlab="PCA axis", ylab="Cumulated variance (%)", cex=1, pch=20, type="h", lwd=2)
}

add.scatter(myInset(), posi="bottomright", inset=c(-0.03,-0.01), ratio=.28, bg=transp("white"))


###################################################
### code chunk number 15: adegenet-dapc.Rnw:420-421
###################################################
scatter(dapc1,1,1, col=myCol, bg="white", scree.da=FALSE, legend=TRUE, solid=.4)


###################################################
### code chunk number 16: adegenet-dapc.Rnw:442-446
###################################################
data(H3N2)
H3N2
pop(H3N2) <- H3N2$other$epid
dapc.flu <- dapc(H3N2, n.pca=30,n.da=10)


###################################################
### code chunk number 17: adegenet-dapc.Rnw:451-453
###################################################
myPal <- colorRampPalette(c("blue","gold","red"))
scatter(dapc.flu, col=transp(myPal(6)), scree.da=FALSE, cell=1.5, cex=2, bg="white",cstar=0)


###################################################
### code chunk number 18: adegenet-dapc.Rnw:457-459
###################################################
set.seed(4)
contrib <- loadingplot(dapc.flu$var.contr, axis=2, thres=.07, lab.jitter=1)


###################################################
### code chunk number 19: adegenet-dapc.Rnw:466-478
###################################################
temp <- seploc(H3N2)
snp906 <- truenames(temp[["906"]])$tab
snp399 <- truenames(temp[["399"]])$tab
freq906 <- apply(snp906, 2, function(e) tapply(e, pop(H3N2), mean, na.rm=TRUE))
freq399 <- apply(snp399, 2, function(e) tapply(e, pop(H3N2), mean, na.rm=TRUE))
freq906
freq399
par(mfrow=c(1,2), mar=c(5.1,4.1,4.1,.1),las=3)
matplot(freq906, pch=c("a","c"), type="b",xlab="year",ylab="allele frequency", xaxt="n", cex=1.5, main="SNP # 906")
axis(side=1, at=1:6, lab=2001:2006)
matplot(freq399, pch=c("c","t"), type="b", xlab="year",ylab="allele frequency", xaxt="n", cex=1.5, main="SNP # 399")
axis(side=1, at=1:6, lab=2001:2006)


###################################################
### code chunk number 20: adegenet-dapc.Rnw:504-507
###################################################
class(dapc1$posterior)
dim(dapc1$posterior)
round(head(dapc1$posterior),3)


###################################################
### code chunk number 21: adegenet-dapc.Rnw:511-512
###################################################
summary(dapc1)


###################################################
### code chunk number 22: adegenet-dapc.Rnw:521-522
###################################################
assignplot(dapc1, subset=1:50)


###################################################
### code chunk number 23: adegenet-dapc.Rnw:539-540
###################################################
compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:6), lab="", ncol=1, xlab="individuals")


###################################################
### code chunk number 24: adegenet-dapc.Rnw:544-545
###################################################
compoplot(dapc1, subset=1:50, posi="bottomright", txt.leg=paste("Cluster", 1:6), lab="", ncol=2, xlab="individuals")


###################################################
### code chunk number 25: adegenet-dapc.Rnw:551-554
###################################################
temp <- which(apply(dapc1$posterior,1, function(e) all(e<0.9)))
temp
compoplot(dapc1, subset=temp, posi="bottomright", txt.leg=paste("Cluster", 1:6),  ncol=2)


###################################################
### code chunk number 26: adegenet-dapc.Rnw:586-589
###################################################
data(microbov)
microbov
temp <- summary(dapc(microbov, n.da=100, n.pca=3))$assign.per.pop*100


###################################################
### code chunk number 27: adegenet-dapc.Rnw:591-593
###################################################
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual breed", horiz=TRUE, las=1)


###################################################
### code chunk number 28: adegenet-dapc.Rnw:601-602
###################################################
temp <- summary(dapc(microbov, n.da=100, n.pca=300))$assign.per.pop*100


###################################################
### code chunk number 29: adegenet-dapc.Rnw:604-606
###################################################
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual breed", horiz=TRUE, las=1)


###################################################
### code chunk number 30: adegenet-dapc.Rnw:614-617
###################################################
x <- microbov
pop(x) <- sample(pop(x))
temp <- summary(dapc(x, n.da=100, n.pca=300))$assign.per.pop*100


###################################################
### code chunk number 31: adegenet-dapc.Rnw:619-621
###################################################
par(mar=c(4.5,7.5,1,1))
barplot(temp, xlab="% of reassignment to actual breed", horiz=TRUE, las=1)


###################################################
### code chunk number 32: adegenet-dapc.Rnw:641-647
###################################################
dapc2 <- dapc(microbov, n.da=100, n.pca=10)
temp <- a.score(dapc2)
names(temp)
temp$tab[1:5,1:5]
temp$pop.score
temp$mean


###################################################
### code chunk number 33: adegenet-dapc.Rnw:651-652
###################################################
dapc2 <- dapc(microbov, n.da=100, n.pca=50)


###################################################
### code chunk number 34: adegenet-dapc.Rnw:654-655 (eval = FALSE)
###################################################
## temp <- optim.a.score(dapc2)


###################################################
### code chunk number 35: adegenet-dapc.Rnw:671-673
###################################################
dapc3 <- dapc(microbov, n.da=100, n.pca=20)
myCol <- rainbow(15)


###################################################
### code chunk number 36: adegenet-dapc.Rnw:675-677
###################################################
par(mar=c(5.1,4.1,1.1,1.1), xpd=TRUE)
compoplot(dapc3, lab="", posi=list(x=12,y=-.01), cleg=.7)


###################################################
### code chunk number 37: adegenet-dapc.Rnw:682-687
###################################################
temp <- which(apply(dapc3$posterior,1, function(e) all(e<0.5)))
temp
lab <- pop(microbov)
par(mar=c(8,4,5,1), xpd=TRUE)
compoplot(dapc3, subset=temp, cleg=.6, posi=list(x=0,y=1.2),lab=lab)


###################################################
### code chunk number 38: adegenet-dapc.Rnw:729-736
###################################################
data(microbov)
set.seed(2)
kept.id <- unlist(tapply(1:nInd(microbov), pop(microbov), function(e) sample(e, 20,replace=FALSE)))
x <- microbov[kept.id]
x.sup <- microbov[-kept.id]
nInd(x)
nInd(x.sup)


###################################################
### code chunk number 39: adegenet-dapc.Rnw:744-750
###################################################
dapc4 <- dapc(x,n.pca=20,n.da=15)
pred.sup <- predict.dapc(dapc4, newdata=x.sup)
names(pred.sup)
head(pred.sup$assign)
pred.sup$ind.scores[1:5,1:3]
round(pred.sup$posterior[1:5, 1:5],3)


###################################################
### code chunk number 40: adegenet-dapc.Rnw:759-767
###################################################
col <- rainbow(length(levels(pop(x))))
col.points <- transp(col[as.integer(pop(x))],.2)
scatter(dapc4, col=col, bg="white", scree.da=0, pch="", cstar=0, clab=0, xlim=c(-10,10), legend=TRUE)
par(xpd=TRUE)
points(dapc4$ind.coord[,1], dapc4$ind.coord[,2], pch=20, col=col.points, cex=5)
col.sup <- col[as.integer(pop(x.sup))]
points(pred.sup$ind.scores[,1], pred.sup$ind.scores[,2], pch=15, col=transp(col.sup,.7), cex=2)
add.scatter.eig(dapc4$eig,15,1,2, posi="bottomright", inset=.02)


###################################################
### code chunk number 41: adegenet-dapc.Rnw:773-774
###################################################
mean(as.character(pred.sup$assign)==as.character(pop(x.sup)))


###################################################
### code chunk number 42: adegenet-dapc.Rnw:780-781
###################################################
table.value(table(pred.sup$assign, pop(x.sup)), col.lab=levels(pop(x.sup)))


