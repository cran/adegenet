\encoding{UTF-8}
\docType{methods}
\name{Auxiliary functions}
\alias{adegenetWeb}
\alias{adegenetTutorial}
\alias{checkType}
\alias{.rmspaces}
\alias{.genlab}
\alias{.readExt}
\title{ Utilities functions for adegenet}
\description{
  These functions are mostly used internally in adegenet.
  The notable exceptions are \code{adegenetWeb} which opens the adegenet
  website in the default navigator, and \code{adegenetTutorial} which
  opens online tutorials for adegenet.\cr

  The other functions are:\cr
  - \code{checkType}: checks the type of markers being used in a
  function and issues an error if appropriate.\cr
  - \code{.rmspaces}: remove peripheric spaces in a character string. \cr
  - \code{.genlab}: generate labels in a correct alphanumeric ordering. \cr
  - \code{.readExt}: read the extension of a given file. \cr
 
}
\usage{
adegenetWeb()
adegenetTutorial(which=c("general","spca"))
.genlab(base, n)
}
\arguments{
  \item{which}{a character string being "general" or "spca", indicating
    which tutorial should be opened.}
  \item{base}{a character string forming the base of the labels}
  \item{n}{the number of labels to generate}
}
\value{
  For \code{.genlab}, a character vector of size "n".
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## this opens the adegenet website
adegenetWeb()

## this opens the general tutorial for adegenet
adegenetTutorial()

## this opens the sPCA tutorial for adegenet
adegenetTutorial("spca")
}

.genlab("Locus-",11)
}
\keyword{manip}