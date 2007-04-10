\name{import}
\alias{import2genind}
\alias{genetix2genind}
\alias{genepop2genind}
\alias{fstat2genind}
\title{ Conversion function for adegenet}
\description{ The function \code{import2genind} detects the extension of
  the file given in argument and seeks for an appropriate import
  function to create a \code{genind} object.\cr
  Current functions are :\cr
  - \code{genetix2genind} for GENETIX files (.gtx). Note that this
  function is called by the others. \cr
  - \code{genepop2genind} for Genepop files (.gen) \cr
  - \code{fstat2genind} for Fstat files \code{.dat} \cr
}
\usage{
import2genind(file,missing=NA,quiet=FALSE)
genetix2genind(file=NULL,X=NULL,pop=NULL,missing=NA,quiet=FALSE)
genepop2genind(file,missing=NA,quiet=FALSE)
fstat2genind(file,missing=NA,quiet=FALSE)
}
\arguments{
  \item{file}{ a character string giving the path to the file to
    convert, with the appropriate extension.}
  \item{missing}{can be NA, 0 or "mean". See details section.}
  \item{quiet}{ logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
  \item{X}{if file is not provided, \code{genetix2genind} can be used on
  a data frame with genotypes in GENETIX format (e.g. "080082" for an
  heterozygote with alleles 80 and 82); individuals are in rows, loci are
  in columns. Missing values are coded as "000000".}
  \item{pop}{an optional factor giving the population of each genotype in 'x'. }}
\details{There are 3 treatments for missing values: \cr
  - NA: kept as NA. \cr
  - 0: missing values are considered as zero. Recommended for a PCA on
  compositionnal data. \cr
  - "mean": missing values are given the mean frequency of the
  corresponding allele. Recommended for a centred PCA.\cr

  Beware: same data in different formats are not expected to produce the
  exactly the same \code{genind} objects.\cr
  For instance, conversions made by GENETIX to Fstat may change the
  the sorting of the genotypes; GENETIX stores individual names whereas
  Fstat does not; Genepop chooses a sample's name from the name of its last
  genotype; etc.
}
\value{an object of the class \code{genind}}
\references{
  Belkhir K., Borsa P., Chikhi L., Raufaste N. & Bonhomme F. (1996-2004)
  GENETIX 4.05, logiciel sous Windows TM pour la génétique des
  populations. Laboratoire Génome, Populations, Interactions, CNRS UMR
  5000, Université de Montpellier II, Montpellier (France). \cr
  
  Raymond M. & Rousset F, (1995). GENEPOP (version 1.2): population
  genetics software for exact tests and ecumenicism. \emph{J. Heredity},
  \bold{86}:248-249 \cr
  
  Fstat (version 2.9.3). Software by Jerome
  Goudet. http://www2.unil.ch/popgen/softwares/fstat.htm\cr

  Excoffier L. & Heckel G.(2006) Computer programs for population genetics data analysis: a survival guide \emph{Nature}, \bold{7}: 745-758
}
\seealso{\code{\link[pkg:hierfstat]{read.fstat.data}}}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
genetix2genind(system.file("files/nancycats.gtx",package="adegenet"))

fstat2genind(system.file("files/nancycats.dat",package="adegenet"))

genepop2genind(system.file("files/nancycats.gen",package="adegenet"))

import2genind(system.file("files/nancycats.gtx",
package="adegenet"))

if(require(hierfstat)){
obj <- fstat2genind(system.file("data/diploid.dat",package="hierfstat"))
obj
}
}
\keyword{manip}