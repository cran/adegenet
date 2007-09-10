\encoding{UTF-8}
\name{genind}
\alias{genind}
\alias{is.genind}
\alias{as.genind}
\alias{print.genind}
\alias{summary.genind}
\alias{genind2genpop}
\alias{.rmspaces}
\alias{.is.gen}
\title{adegenet class for individual genotypes}
\description{The objects of class \code{genind} contain individual genotypes.\cr
  It consists in a list with several components (see value section).\cr
  The function \code{genind2genpop} converts individuals genotypes of
  known population into a \code{genpop} object.\cr
  The \code{summary} of a \code{genind} object invisibly returns a list of components
  (see value section). 
  The function \code{as.genind} is called by import functions (see \code{\link{import2genind}}).
}
\usage{
is.genind(x)
as.genind(tab=NULL,pop=NULL,prevcall=NULL)
\method{print}{genind}(x, \dots)
\method{summary}{genind}(object, \dots)
genind2genpop(x,pop=NULL,missing=NA,quiet=FALSE)
}
\arguments{
  \item{x}{an object of class \code{genind}.}
  \item{tab}{a individuals x alleles matrix of genotypes coded as
    allelic frequencies.}
  \item{pop}{a factor giving the population of each genotype in 'x'. If
    note provided, seeked in x\$pop, but if given, the argument prevails on x\$pop.}
  \item{prevcall}{call of an object, for internal use.}
  \item{...}{other -unused- arguments}
  \item{object}{an object of class \code{genind}.}
  \item{missing}{can be NA, 0, or "replace". See details for more information.}
  \item{quiet}{logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
}
\details{The values of the 'missing' argument in \code{genind2genpop}
  have the following effects:\cr
  - NA: if all genotypes of a population for a given allele are
  missing, count value will be NA\cr
  - 0: if all genotypes of a population for a given allele are
  missing, count value will be 0\cr
  - "replace": when an allele is not typed in a population, it is
  assigned an allele count so that the allelic frequency in this
  populations is the same as the frequency in the whole dataset.\cr
  If allele 'j' of locus 'k' in pop 'i' is missing, the
  count value is number 'x' so that the frequency 'x/s' ('s' being the
  number of observations in 'k' ) equals the frequency 'f' computed on the whole data (i.e. considering all pop as one)\cr
  Then x verifies:\cr
  \deqn{x/s = f(1-f) => x=f(1-f)s}

  Note on 'pop': the factor is sorted in the order of appearance of each modality.
}
\value{
  \item{tab}{matrix of genotypes -in rows- for all alleles -in
    columns-. Values are frequency: '0' if the genotype does not have
    the corresponding allele, '1' for an homozygote and 0.5 for an heterozygte.Rows and columns are given generic names.}
  \item{ind.names}{character vector containing the real names of the
    individuals. Note that as Fstat does not store these names, objects
    converted from .dat files will contain empty \code{ind.names}.}
  \item{loc.names}{character vector containing the real names of the loci}
  \item{loc.nall}{integer vector giving the number of alleles per locus}
  \item{loc.fac}{locus factor for the columns of \code{tab}}
  \item{all.names}{list having one component per locus, each containing a character vector of alleles names}
  \item{call}{the matched call}
  \item{pop}{(optional) factor giving the population of each individual}
  \item{pop.names}{(optional) vector giving the real names of the
    populations}
  \item{N}{(summary) total number of genotypes.}
  \item{pop.eff}{(summary) populations sample size.}
  \item{loc.nall}{(summary) number of alleles per locus.}
  \item{pop.nall}{(summary) number of alleles per population.}
  \item{NA.perc}{(summary) percentage of - appearing - missing data.}
  \item{Hobs}{(summary) observed heterozygosity.}
  \item{Hexp}{(summary) expected heterozygosity.}  
}
\references{}
\seealso{\code{\link{genpop}}, \code{\link{import2genind}}, \code{\link{genetix2genind}}, 
\code{\link{genepop2genind}}, \code{\link{fstat2genind}}
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
obj <- genetix2genind(system.file("files/nancycats.gtx",package="adegenet"),missing="mean")
is.genind(obj)
summary(obj)
obj

# test inter-colonies structuration
if(require(hierfstat)){
gtest <- gstat.randtest(obj,nsim=99)
gtest
plot(gtest)
}

# perform an inter-class PCA
if(require(ade4)){
pca1 <- dudi.pca(obj$tab,scannf=FALSE,scale=FALSE)
pcabet1 <- between(pca1,obj$pop,scannf=FALSE)
pcabet1

s.class(pcabet1$ls,obj$pop,sub="Inter-class PCA",possub="topleft",csub=2)
add.scatter.eig(pcabet1$eig,2,xax=1,yax=2)
}
}
\keyword{manip}
\keyword{multivariate}

