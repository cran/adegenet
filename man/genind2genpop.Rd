\encoding{UTF-8}
\name{genind2genpop}
\alias{genind2genpop}
\title{Conversion from a genind to a genpop object}
\description{The function \code{genind2genpop} converts genotypes data
  (genind) into alleles counts per population (genpop).}
\usage{
genind2genpop(x,pop=NULL,missing=c("NA","0","chi2"),quiet=FALSE)
}
\arguments{
  \item{x}{an object of class \code{genind}.}
  \item{pop}{a factor giving the population of each genotype in 'x'. If
    note provided, seeked in x@pop, but if given, the argument
    prevails on x@pop.}
  \item{missing}{can be "NA", "0", or "chi2". See details for more information.}
  \item{quiet}{logical stating whether a conversion message must be
    printed (TRUE,default) or not (FALSE).}
    
}
\details{The values of the 'missing' argument in \code{genind2genpop}
  have the following effects:\cr
  - "NA": if all genotypes of a population for a given allele are
  missing, count value will be NA\cr
  
  - "0": if all genotypes of a population for a given allele are
  missing, count value will be 0\cr
  
  - "chi2": if all genotypes of a population for a given allele are
  missing, count value will be that of a theoretical count in of a
  Chi-squared test. This is obtained by the product of the margins
  sums divided by the total number of alleles.
}
\value{
  A genpop object. The component @other in 'x' is passed to the created
  genpop object.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\seealso{
  \linkS4class{genind}, \linkS4class{genpop}, \code{\link{na.replace}}
}
\examples{
  data(nancycats)
  nancycats
  catpop <- genind2genpop(nancycats)
  catpop
  summary(catpop)
}
\keyword{classes}
\keyword{manip}
\keyword{multivariate}
