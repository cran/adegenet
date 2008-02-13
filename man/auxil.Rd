\encoding{UTF-8}
\docType{methods}
\name{Auxiliary functions}
\alias{adegenetWeb}
\alias{.rmspaces}
\alias{.genlab}
\alias{.readExt}
\title{ Utilities functions for adegenet}
\description{
  These functions are mostly used internally in adegenet.
  The notable exception is \code{adegenetWeb}, which opens the adegenet
  website in the default navigator.\cr

  The other functions are:\cr
  - \code{.rmspaces}: remove peripheric spaces in a character string \cr
  - \code{.genlab}: generate labels in a correct alphanumeric ordering \cr
  - \code{.readExt}: read the extension of a given file \cr
 
}
\usage{
adegenetWeb()
.genlab(base, n)
}
\arguments{
    \item{base}{a character string forming the base of the labels}
    \item{n}{the number of labels to generate}
}
\value{
  For \code{.genlab}, a character vector of size "n".
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
\dontrun{
## this opens the adegenet website
adegenetWeb()
}

.genlab("Locus-",11)
}
\keyword{manip}