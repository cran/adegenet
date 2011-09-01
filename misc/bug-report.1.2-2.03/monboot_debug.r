## monboot.r
## An R program designed to take geographic coordinates and a genetic distance matrices produced by a companion C++ program (monboot.cxx)
## and compute a Monmonier's barrier.  
## In addition, the program is designed to take many bootstrapped matrices for analysis.
## The goal of the bootstrap is to find the highly likely region where there is a major genetic break on the landscape.

## Algorithm
#	Compute Monmonier barrier on original data.  Output that info.
#		Read in original coordinates.  Store.
#		Read in genetic distance matrix.  Store.
#		Compute network.
#		compute barrier.
#		Output barrier.
#	Perform bootstrap analysis. (loop)
#		Input filenames are predetermined based on companion C++ output.
#		Read in coordinates.  Store.
#		Read in genetic distance matrix.  Store.
#		Compute network.
#		compute barrier.
#		Output barrier.

#############################

#### User defined variables (to be changed to appropriate when used)

pathname <- "monmonierwork/EV/EVmonboot/"  #base path for directory that all input and output data will be in
#numboot <- 1000  #number of bootstrap replicates to run
b <- 7

## Input files
orig_coord_filename <- paste(pathname, "orig_coord.input", sep = "")
boot_coord_base_filename <- paste(pathname, "boot_coord.", sep = "")
orig_dist_filename <- paste(pathname, "orig_dist.input", sep ="")
boot_dist_base_filename <- paste(pathname, "boot_dist.", sep = "")

## Output files
orig_mon_filename <- paste(pathname, "orig_mon.out", sep = "")
boot_mon_filename <- paste(pathname, "boot_mon.out", sep = "")

#### Necessary libraries

library(spdep)
library(ade4)
library(adegenet)

#### Main program

## 1. Computing Monmonier barrier on original data.

# Read in coordinates and genetic matrix
orig_coords <- read.table(orig_coord_filename)  
numpops <- length(orig_coords[,1])
orig_Da <- read.table(orig_dist_filename)
orig_Da <- as.dist(orig_Da)

# Compute the barrier
orig_cn <- chooseCN(orig_coords, ask=FALSE, type=1)
orig_mon <- optimize.monmonier(coords, Da, cn, threshold=0)
orig_mon <- monmonier(orig_coords, orig_Da, orig_cn, threshold=0, nrun=1, skip.local.diff = rep(1, 1)) #threshold is arbitrarily set at 0
plot(orig_mon)

# Output barrier info -- for now only outputting data from run1
write.table(orig_mon$run1$dir1$path, file=orig_mon_filename, sep="\t")
write.table(orig_mon$run1$dir2$path, file=orig_mon_filename, sep="\t", append = TRUE, col.names = FALSE)

## 2. Performing the bootstrap analyses

#clear the output files
write("",file=boot_mon_filename)

#for (b in 1:numboot)
#{
	# Read in coordinates and genetic matrix
	coords <- read.table(paste(boot_coord_base_filename, b, ".input", sep=""))
	numpops <- length(coords[,1])
	Da <- read.table(paste(boot_dist_base_filename, b, ".input", sep=""))
	Da <- as.dist(Da)

	# Calc monmonier barrier
	cn <- chooseCN(coords, ask=FALSE, type=1)
	mon <- monmonier(coords, Da, cn, threshold=0, nrun=1, skip.local.diff = rep(1, 1)) #threshold is arbitrarily set at 0
	mon <- optimize.monmonier(coords, Da, cn, threshold=0)
#	plot(mon)

	# Output info
	write.table(cbind(rep(b), mon$run1$dir1$path), file=boot_mon_filename, sep="\t", append = TRUE, col.names = FALSE)
	write.table(cbind(rep(b), mon$run1$dir2$path), file=boot_mon_filename, sep="\t", append = TRUE, col.names = FALSE)
#	write("\n",file=boot_mon_filename, sep="\t", append = TRUE)
#}
