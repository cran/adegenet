\encoding{UTF-8}
\name{genind constructor}
\alias{genind-methods}
\alias{genind}
\alias{as.genind}
\alias{is.genind}
\title{genind constructor}
\description{Constructor for \linkS4class{genind} objects.\cr
  The function \code{genind} creates a \linkS4class{genind} object from
  a matrix of allelic frequency where genotypes are in rows and alleles
  in columns. This table must have correct names for rows and
  columns.\cr

  The function \code{as.genind} is an alias for \code{genind} function.\cr
 
  \code{is.genind} tests if an object is a valid genind object.\cr

  Note: to get the manpage about \linkS4class{genind}, please type
  'class ? genind'.
}

\usage{
genind(tab,pop=NULL,prevcall=NULL)
as.genind(tab,pop=NULL,prevcall=NULL)
is.genind(x)
}

\arguments{
  \item{tab}{a individuals x alleles matrix of genotypes coded as
    allelic frequencies, i.e. like in a genind object}
  \item{pop}{a factor giving the population of each genotype in 'x'}
  \item{prevcall}{call of an object}
  \item{x}{an object}
}
\value{
  For \code{genind} and \code{as.genind}, a genind object. For \code{is.genind}, a logical.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr}}
\seealso{
  \code{"\linkS4class{genind}"}, and \code{\link{import2genind}} for
  importing from various types of file.
}
\examples{
data(nancycats)
nancycats@loc.names

# isolate one marker, fca23
obj <- seploc(nancycats)$"fca23"
obj
}
\keyword{manip}
