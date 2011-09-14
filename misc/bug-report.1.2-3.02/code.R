## PROBLEM 1: ISSUE WITH NON-NUMERIC ALLELES

## This simple example already triggers a problem
df=matrix(c('a/b','a/a','a/b','x/x','NA','x/y'),nrow=2)
colnames(df)=paste("locus",1:3,sep=".")
rownames(df)=1:2
df # looks ok

toto=df2genind(df, sep="/", ploidy=2)
toto

## here is the issue
propShared(toto)
genind2df(toto,sep="/")





## PROBLEM 2:...
dat <- read.csv("example.csv")
rownames(dat) <- dat[,1]
dat <- dat[,-1]
x <- df2genind(dat, sep="\\+")
propShared(x)
