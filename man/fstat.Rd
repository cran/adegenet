\encoding{UTF-8}
\name{F statistics}
\alias{fstat}
\alias{Fst}
\alias{FST}
\alias{fst}
\alias{pairwise.fst}
\title{F statistics for genind objects}
\description{
  The function \code{fstat} computes a global Fst, while
  \code{pairwise.fst} computes Nei's pairwise Fst between all pairs of
  populations. Both functions are designed for \linkS4class{genind}
  objects.

  \code{fstat} is wrapper for \code{varcomp.glob} from package
  \code{hierfstat} for \linkS4class{genind} objects. It computes F statistics (Fst, Fis, Fit)
  given a set of genotypes and a grouping factor.

  \code{pairwise.fst} is an implementation of Nei's Fst in which
  heretozygosities are weighted by group sizes (see details).
}
\usage{
fstat(x, pop=NULL, fstonly=FALSE)

pairwise.fst(x, pop=NULL, res.type=c("dist","matrix"), truenames=TRUE)
}
\arguments{
  \item{x}{an object of class \linkS4class{genind}.}
  \item{pop}{a factor giving the 'population' of each individual. If NULL,
    pop is seeked from \code{pop(x)}. Note that the term population refers in
    fact to any grouping of individuals'.}
  \item{fstonly}{a logical stating whether only the Fst value should be
    returned (TRUE) instead of all F statistics (FALSE, default).}
  \item{res.type}{the type of result to be returned: a \code{dist} object, or a
    symmetric matrix}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
}
\value{
  A vector, a matrix, or a dist object containing F statistics.}
\seealso{\code{\link{Hs}}, \code{\link[hierfstat]{varcomp.glob}},
  \code{\link{gstat.randtest}}
}
\references{
  Nei, M. (1973) Analysis of gene diversity in subdivided
  populations. Proc Natl Acad Sci USA, 70: 3321-3323
}
\details{
  Let \eqn{A} and \eqn{B} be two populations of population sizes \eqn{n_A} and
  \eqn{n_B}, with expected heterozygosity (averaged over loci) \eqn{Hs(A)} and \eqn{Hs(B)},
  respectively. We denote \eqn{Ht} the expected heterozygosity of a population
  pooling \eqn{A} and \eqn{B}. Then, the pairwise \eqn{Fst} between \eqn{A} and \eqn{B} is computed
  as:\cr

 \eqn{  Fst(A,B) = \frac{(Ht - (n_A Hs(A) + n_B Hs(B))/(n_A + n_B) )}{Ht}} \cr
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
if(require(hierfstat)){
data(nancycats)

## Fst, Fis, Fit
fstat(nancycats)

## pairwise Fst
mat.fst <- pairwise.fst(nancycats, res.type="matrix")
mat.fst
}
}
\keyword{multivariate}