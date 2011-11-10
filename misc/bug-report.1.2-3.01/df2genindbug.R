df=matrix(c('a/b','a/a','a/b','x/x','NA','x/y'),nrow=2)
colnames(df)=paste("locus",1:3,sep=".")
rownames(df)=1:2
df # looks ok

toto=df2genind(df, sep="/", ploidy=2)
toto
toto@tab # only one locus !?!
genind2df(toto,sep="/")
