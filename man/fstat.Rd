\encoding{UTF-8}
\name{fstat}
\alias{fstat}
\title{F statistics for genind objects}
\description{This function is a wrapper of \code{varcomp.glob} for
  \linkS4class{genind} objects. It computes F statistics (Fst, Fis, Fit)
  given a set of genotypes and a grouping factor.
}
\usage{
fstat(x, pop=NULL, fstonly=FALSE)
}
\arguments{
  \item{x}{an object of class \linkS4class{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
  pop is seeked from x\@pop. Note that the term population refers in
  fact to any grouping of individuals'.}
\item{fstonly}{a logical stating whether only the Fst value should be
  returned (TRUE) instead of all F statistics (FALSE, default).}
  }
\value{A matrix of F statistics.}
\seealso{\code{\link[hierfstat]{varcomp.glob}},
  \code{\link{gstat.randtest}}
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
if(require(hierfstat)){
data(nancycats)
fstat(nancycats)
}
}
\keyword{multivariate}