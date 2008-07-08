\name{scaleGen-methods}
\docType{methods}
\alias{scaleGen}
\alias{scaleGen-methods}
\alias{scaleGen,genind-method}
\alias{scaleGen,genpop-method}
\title{ Compute scaled allele frequencies }
\description{
  The generic function \code{scaleGen} is an analogue to the
  \code{scale} function, but is designed with further arguments giving
  scaling options.\cr

  Methods are defined for \linkS4class{genind} and \linkS4class{genpop}
  objects.
  Both return data.frames of scaled allele frequencies.
}
\usage{
\S4method{scaleGen}{genind}(x, center=TRUE, scale=TRUE, method=c("sigma", "binom"), missing=c("NA","0","mean"),truenames=TRUE)
\S4method{scaleGen}{genpop}(x, center=TRUE, scale=TRUE, method=c("sigma", "binom"), missing=c("NA","0","mean"),truenames=TRUE)
}
\arguments{
  \item{x}{a \linkS4class{genind} and \linkS4class{genpop} object}
  \item{center}{a logical stating whether alleles frequencies should be
    centred to mean zero (default to TRUE). Alternatively, a vector of
    numeric values, one per allele, can be supplied: these values will be
    substracted from the allele frequencies.}
  \item{scale}{a logical stating whether alleles frequencies should be
    scaled (default to TRUE). Alternatively, a vector of
    numeric values, one per allele, can be supplied: these values will be
    substracted from the allele frequencies.}
  \item{method}{a character indicating the method to be used. See details.}
  \item{truenames}{a logical indicating whether true labels (as opposed
    to generic labels) should be used to name the output.}
 \item{missing}{a character giving the treatment for missing values. Can be "NA", "0" or "mean"}
 }
 \value{
   A \linkS4class{genind} and \linkS4class{genpop} object.
}
\details{
  The argument \code{method} is used as follows:\cr

  - \code{sigma}: scaling is made using the usual standard deviation\cr

  - \code{binom}: scaling is made using the theoretical variance of the
  allele frequency. This can be used to avoid that frequencies close to
  0.5 have a stronger variance that those close to 0 or 1.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{

}
\keyword{methods}
\keyword{manip}
