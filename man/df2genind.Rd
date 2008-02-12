\encoding{UTF-8}
\name{df2genind}
\alias{df2genind}
\alias{genind2df}
\title{ Convert a data.frame of genotypes to a genind object, and conversely.}
\description{
  The function \code{df2genind} converts a data.frame (or a matrix) into
  a \linkS4class{genind} object. The data.frame must meet the following
  requirements:\cr
  - genotypes are in row (on row per genotype)\cr
  - markers are in columns\cr
  - each element is a string of (2, 4 or 6) characters coding both alleles without
  separation. Missing values are coded by NA or zeros. Uncomplete
  strings are filled with "0" at the begining. For instance, "912" in a
  4-characters coding scheme is interpreted as "0912", thus as an
  heterozygote with alleles "09" and "12". Note that this format is the
  one of GENETIX for a 6-characters coding scheme\cr

  The function \code{genind2df} converts a \linkS4class{genind} back to
 such a data.frame.
 }
\usage{
df2genind(X, ncode=NULL, ind.names=NULL, loc.names=NULL, pop=NULL,
 missing=NA)
genind2df(x,pop=NULL, sep="")
}
\arguments{
  \item{X}{a matrix or a data.frame (see decription)}
  \item{ncode}{an optional integer the number of characters used for coding one
    genotype at one locus (can be 2, 4 or 6). If not provided, this is
    determined from data.}
  \item{ind.names}{an optional character vector giving the individuals
    names; if NULL, taken from rownames of X.}
  \item{loc.names}{an optional character vector giving the markers
    names; if NULL, taken from colnames of X.}
  \item{pop}{an optional factor giving the population of each individual.}
  \item{missing}{can be NA, 0 or "mean". See details section.}
  \item{x}{a \linkS4class{genind} object}
  \item{sep}{a character used to separate two alleles}
 }
   
\details{There are 3 treatments for missing values: \cr
  - NA: kept as NA.\cr
  
  - 0: allelic frequencies are set to 0  on all alleles of the concerned
  locus. Recommended for a PCA on compositionnal data.\cr
  
  - "mean": missing values are replaced by the mean frequency of the
  corresponding allele, computed on the whole set of
  individuals. Recommended for a centred PCA.\cr
}
\value{an object of the class \linkS4class{genind} for \code{df2genind};
  a matrix of biallelic genotypes for \code{genind2df}}
\seealso{\code{\link{import2genind}}, \code{\link{read.genetix}},
  \code{\link{read.fstat}}, \code{\link{read.structure}}}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
df <- data.frame(loc1=c("11","11","12","32"),
loc2=c(NA,"34","55","15"),loc3=c("22","22","21","22"))
row.names(df) <- .genlab("geontype",4)
df

obj <- df2genind(df)
obj
truenames(obj)
genind2df(obj)
}
\keyword{manip}
