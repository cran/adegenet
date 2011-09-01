##
## makefreq doesn't work after seploc
##
data(nancycats)
x <- seploc(genind2genpop(nancycats))$fca37
makefreq(x)
