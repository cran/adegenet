\name{haploGen}
\alias{haploGen}
\alias{print.haploGen}
\alias{[.haploGen}
\alias{labels.haploGen}
\alias{as.POSIXct.haploGen}
\alias{seqTrack.haploGen}
\alias{haploGen-class}
% \alias{seqTrackG.haploGen}
% \alias{optimize.seqTrack.haploGen}
\alias{as.seqTrack.haploGen}
\alias{plotHaploGen}
\alias{sample.haploGen}
%\alias{as,haploGen,graphNEL-method}
%\alias{coerce,haploGen,graphNEL-method}
%\alias{as,haploGen,graphNEL-method}
\title{Simulation of genealogies of haplotypes}
\description{
  The function \code{haploGen} implements simulations of genealogies of
  haplotypes. This forward-time, individual-based simulation tool allow
  haplotypes to replicate and mutate according to specified parameters,
  and keeps track of entire genealogies.

  Simulations can be spatially explicit or not (see \code{geo.sim}
  argument). In the first case, haplotypes are assigned to locations on
  a regular grip. New haplotypes disperse from their ancestor's location
  according to a random Poisson diffusion, or alternatively according to
  a pre-specified migration scheme. This tool does not allow for
  simulating selection or linkage disequilibrium.

  Produced objects are lists with the class \code{haploGen}; see 'value'
  section for more information on this class. Other functions are
  available to print, plot, subset, sample or convert \code{haploGen}
  objects. A seqTrack method is also provided for analysing
  \code{haploGen} objects.
}
\usage{
haploGen(seq.length=10000, mu=0.0001, t.max=20,
              gen.time=function(){round(rnorm(1,5,1))},
              repro=function(){round(rnorm(1,2,1))}, max.nb.haplo=1e3,
              geo.sim=TRUE, grid.size=5, lambda.xy=0.5,
              mat.connect=NULL, ini.n=1, ini.xy=NULL)
\method{print}{haploGen}(x, \dots)
\method{[}{haploGen}(x, i, j, drop=FALSE)
\method{labels}{haploGen}(object, \dots)
\method{as.POSIXct}{haploGen}(x, tz="", origin=as.POSIXct("2000/01/01"), \dots)
\method{seqTrack}{haploGen}(x, best=c("min","max"), prox.mat=NULL, \dots)
as.seqTrack.haploGen(x)
plotHaploGen(x, annot=FALSE, date.range=NULL, col=NULL, bg="grey", add=FALSE, \dots)
sample.haploGen(x, n)
%\S4method{coerce}{haploGen,graphNEL}(from, to, strict=TRUE)
}
\arguments{
  \item{seq.length}{an integer indicating the length of the simulated
    haplotypes, in number of nucleotides.}
  \item{mu}{the mutation rate, in number of mutation per site and per
    time unit. Can be a (fixed) number or a function returning a number
    (then called for each replication event).}
  \item{t.max}{an integer indicating the maximum number of time units to
    run the simulation for.}
  \item{gen.time}{an integer indicating the generation time, in number
    of time units. Can be a (fixed) number or a function returning a number
    (then called for each reproduction event).}
  \item{repro}{an integer indicating the number of descendents per
    haplotype. Can be a (fixed) number or a function returning a number
    (then called for each reproduction event).}
  \item{max.nb.haplo}{an integer indicating the maximum number of
    haplotypes handled at any time of the simulation, used to control the
    size of the produced object. Larger number will lead to slower
    simulations. If this number is exceeded, the genealogy is prunded to as
    to keep this number of haplotypes.}
  \item{geo.sim}{a logical stating whether simulations should be
    spatially explicit (TRUE, default) or not (FALSE). Spatially-explicit
    simulations are slightly slower than their non-spatial counterpart.}
  \item{grid.size}{the size of the square grid of possible locations for
    spatial simulations. The total number of locations will be this number
    squared.}
  \item{lambda.xy}{the parameter of the Poisson distribution used to
    determine dispersion in x and y axes.}
  \item{mat.connect}{a matrix of connectivity describing migration
    amongts all pairs of locations. \code{mat.connect[i,j]} indicates the
    probability, being in 'i', to migrate to 'j'. The rows of this matrix
    thus sum to 1. It has as many rows and columns as there are locations,
    with row 'i' / column 'j' corresponding to locations number 'i' and 'j'.
    Locations are numbered as in a matrix in which rows and columns are
    respectively x and y coordinates. For instance, in a 5x5 grid, locations
    are numbered as in \code{matrix(1:25,5,5)}.
  }
  \item{ini.n}{an integer specifying the number of (identical)
    haplotypes to initiate the simulation}
  \item{ini.xy}{a vector of two integers giving the x/y coordinates of the initial haplotype.}
  \item{x,object}{\code{haploGen} objects.}
  \item{i,j, drop}{\code{i} is a vector used for subsetting the object. For
    instance, \code{i=1:3} will retain only the first three haplotypes of the
    genealogy. \code{j} and \code{drop} are only provided for compatibility,
    but not used.}
  \item{best, prox.mat}{arguments to be passed to the
    \code{\link{seqTrack}} function. See documentation of
    \code{\link{seqTrack}} for more information.}
  \item{annot,date.range,col,bg,add}{arguments to be passed to \code{\link{plotSeqTrack}}.}
  \item{n}{an integer indicating the number of haplotypes to be retained
    in the sample}
  % \item{from, to}{arguments of the conversion function, for converting a
  %   \code{haploGen} object into a \code{graphNEL-class}.}
  \item{tz, origin}{ aguments to be passed to \code{\link{as.POSIXct}}
    (see ?as.POSIXct)}
  \item{\dots}{further arguments to be passed to other methods}
  % \item{strict}{a logical used for compatibility with \code{as} generic
  %   function, but not used in the conversion. See \code{\link{setAs}} for
  %   more information.}
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk}}
\references{
  Jombart T, Eggo R, Dodd P, Balloux F (2010) Reconstructing disease
  outbreaks from genetic data: a graph approach. Heredity. doi:
  10.1038/hdy.2010.78.
}
\value{
  === haploGen class ===\cr
  \code{haploGen} objects are lists containing the following slots:\cr
  - seq: DNA sequences in the DNAbin matrix format\cr
  - dates: dates of appearance of the haplotypes\cr
  - ances: a vector of integers giving the index of each haplotype's
  ancestor\cr
  - id: a vector of integers giving the index of each haplotype\cr
  - xy: (optional) a matrix of spatial coordinates of haplotypes\cr
  - call: the matched call


  === misc functions ===\cr
  - as.POSIXct: returns a vector of dates with POSIXct format\cr
  - labels: returns the labels of the haplotypes\cr
  - as.seqTrack: returns a seqTrack object. Note that this object is not a
  proper seqTrack analysis, but just a format conversion convenient for
  plotting \code{haploGen} objects.
}
\details{
  === Dependencies with other packages ===\cr
  - ape package is required as it implements efficient handling of DNA
  sequences used in \code{haploGen} objects. To install this package,
  simply type:\cr
  \code{install.packages("ape")}

  - for various purposes including plotting, converting genealogies to
  graphs (\code{graphNEL-class} class) can be useful. This requires the
  packages graph, and possibly Rgraphviz for plotting. These packages are
  not on CRAN, but on Bioconductor. To install them, use:\cr
  source("http://bioconductor.org/biocLite.R")\cr
  biocLite("graph")\cr
  biocLite("Rgraphviz")

  See the respective vignettes for more information on using these packages.


  % === Converting haploGen objects to graphs ===\cr
  % \code{haploGen} objects can be converted to \code{graphNEL-class}
  % objects, which can in turn be plotted and manipulated using classical
  % graph tools. Simply use 'as(x, "graphNEL")' where 'x' is a
  % \code{haploGen} object. This functionality requires the \code{graph}
  % package (see 'details').
}
\examples{
\dontrun{
if(require(ape)){
## PERFORM SIMULATIONS
x <- haploGen(repro=2)
x

## PLOT SPATIAL SPREAD
plotHaploGen(x, bg="white")
title("Spatial dispersion of the haplotypes")

% ## PLOT GENEALOGY
% if(require(graph) & require(Rgraphviz)){
% g=as(x, "graphNEL")
% g
% renderGraph(layoutGraph(g))
% }


## USE SEQTRACK RECONSTRUCTION
x.recons <- seqTrack(x)
mean(x.recons$ances==x$ances, na.rm=TRUE) # proportion of correct reconstructions
}
}
}
