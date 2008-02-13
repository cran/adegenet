\encoding{UTF-8}
\docType{methods}
\name{seploc}
\alias{seploc}
\alias{seploc-methods}
\alias{seploc,ANY-method}
\alias{seploc,genind-method}
\alias{seploc,genpop-method}
\title{ Separate data per locus}
\description{
  The function \code{seploc} splits an object (\linkS4class{genind} or
  \linkS4class{genpop}) by marker, returning a list of objects whose components
  each correspond to a marker.\cr
}
\usage{
\S4method{seploc}{genind}(x,truenames=TRUE,res.type=c("genind","matrix"))
\S4method{seploc}{genpop}(x,truenames=TRUE,res.type=c("genpop","matrix"))
}
\arguments{
    \item{x}{a \linkS4class{genind} or a \linkS4class{genpop} object.}
    \item{truenames}{a logical indicating whether true names should be
      used (TRUE, default) instead of generic labels (FALSE).}
    \item{res.type}{a character indicating the type of returned results,
      a genind or genpop object (default) or a matrix of data
      corresponding to the 'tab' slot.}
}
\value{The function \code{seploc} returns an list of objects of the
  same class as the initial object, or a list of matrices similar to
  x\$tab.\cr
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\seealso{\code{\link{seppop}}, \code{\link{repool}}}
\examples{
data(microbov)

# separate all markers
obj <- seploc(microbov)
names(obj)

obj$INRA5
}
\keyword{manip}