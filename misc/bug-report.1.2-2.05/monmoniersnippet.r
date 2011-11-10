## Input files
boot_coord_filename <- "boot_coord.21.input"
boot_dist_filename <- "boot_dist.21.input"

#### Necessary libraries

library(spdep)
library(ade4)
library(adegenet)
library(gplots)
source("import.R")
source("monmonier.R")

#### Code snippet

coords <- read.table(boot_coord_filename)
Da <- read.table(boot_dist_filename)
Da <- as.dist(Da)

# Calc monmonier barrier
cn <- chooseCN(coords, ask=FALSE, type=1, plot.nb = FALSE)
mon <- monmonier(coords, Da, cn, threshold=NULL, scanthres = 0, nrun=1) #threshold is arbitrarily set at 0
