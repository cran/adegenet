\encoding{UTF-8}
\docType{methods}
\name{seppop}
\alias{seppop}
\alias{seppop-methods}
\alias{seppop,ANY-method}
\alias{seppop,genind-method}
\title{ Separate genotypes per population}
\description{
  The function \code{seppop} splits a \linkS4class{genind} object by
  population, returning a list of objects whose components each
  correspond to a population.\cr

  By default, components of the list are \linkS4class{genind} objects.
  It can also be a matrix of genotypes corresponding to the x\$tab.\cr
}
\usage{
\S4method{seppop}{genind}(x,pop=NULL,truenames=TRUE,res.type=c("genind","matrix"),
  drop=FALSE, treatOther=TRUE)
}
\arguments{
  \item{x}{a \linkS4class{genind} object}
  \item{pop}{a factor giving the population of each genotype in 'x'. If
    not provided, seeked in x\$pop.}
  \item{truenames}{a logical indicating whether true names should be
    used (TRUE, default) instead of generic labels (FALSE); used if
    res.type is "matrix".}
  \item{res.type}{a character indicating the type of returned results,
    a list of \linkS4class{genind} object (default) or a matrix of data
    corresponding to the 'tab' slots.}
  \item{drop}{a logical stating whether alleles that are no longer
    present in a subset of data should be discarded (TRUE) or kept anyway
    (FALSE, default).}
  \item{treatOther}{a logical stating whether elements of the
  \code{@other} slot should be treated as well (TRUE), or not
  (FALSE). See details in accessor documentations (\code{\link{pop}}).}
}
\value{
  According to 'res.type': a list of \linkS4class{genind} object
  (default) or a matrix of data corresponding to the 'tab' slots.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\seealso{\code{\link{seploc}}, \code{\link{repool}}}
\examples{
data(microbov)

obj <- seppop(microbov)
names(obj)

obj$Salers
}
\keyword{manip}