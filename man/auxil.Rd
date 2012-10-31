\encoding{UTF-8}
\docType{methods}
\name{Auxiliary functions}
\alias{adegenetWeb}
%\alias{adegenetTutorial}
\alias{checkType}
\alias{.rmspaces}
\alias{.genlab}
\alias{.readExt}
\alias{corner}
\alias{num2col}
\alias{fac2col}
\alias{transp}
\alias{bluepal}
\alias{redpal}
\alias{greenpal}
\alias{flame}
\alias{seasun}
\alias{lightseasun}
\alias{deepseasun}
\alias{funky}

\title{ Auxiliary functions for adegenet}
\description{
  adegenet implements a number of auxiliary procedures that might be of
  interest for users. These include graphical tools to translate
  variables (numeric or factors) onto a color scale, adding transparency
  to existing colors, pre-defined color palettes, extra functions to
  access documentation, and low-level treatment of character vectors.
  
  These functions are mostly auxiliary procedures used internally in
  adegenet, with the exception of, which opens the
  adegenet website in the default navigator.\cr

  These items include:\cr
  - \code{adegenetWeb}: opens the adegenet website in a web navigator
  - \code{num2col}: translates a numeric vector into colors. \cr
  - \code{fac2col}: translates a numeric vector into colors. \cr
  - \code{transp}: adds transparency to a vector of colors. Note that
  transparent colors are not supported on some graphical devices.\cr
  - \code{corner}: adds text to a corner of a figure. \cr
  - \code{checkType}: checks the type of markers being used in a
  function and issues an error if appropriate.\cr
  - \code{.rmspaces}: remove peripheric spaces in a character string. \cr
  - \code{.genlab}: generate labels in a correct alphanumeric ordering. \cr
  - \code{.readExt}: read the extension of a given file. \cr

  Color palettes include:\cr
  - \code{bluepal}: white->blue\cr
  - \code{redpal}: white->red\cr
  - \code{greenpal}: white->green\cr
  - \code{flame}: gold->red\cr
  - \code{seasun}: blue->gold->red\cr
  - \code{lightseasun}: blue->gold->red (light variant)\cr
  - \code{deepseasun}: blue->gold->red (deep variant)\cr
  - \code{funky}: many colors\cr
}
\usage{
adegenetWeb()
.genlab(base, n)
corner(text, posi="topleft",  inset=0.1, \dots)
num2col(x, col.pal=heat.colors, reverse=FALSE,
               x.min=min(x), x.max=max(x), na.col="green")
fac2col(x, col.pal=funky, na.col="grey", seed=1)
transp(col, alpha=.5)
}
\arguments{
  \item{base}{a character string forming the base of the labels}
  \item{n}{the number of labels to generate}
  \item{text}{a character string to be added to the plot}
  \item{posi}{a character matching any combinations of "top/bottom" and
    "left/right".}
  \item{inset}{a vector of two numeric values (recycled if needed)
    indicating the inset, as a fraction of the plotting region.}
  \item{\dots}{further arguments to be passed to \code{\link{text}}}
  \item{x}{a numeric vector (for \code{num2col}) or a vector converted
    to a factor (for \code{fac2col}).}
  \item{col.pal}{a function generating colors according to a given palette.}
  \item{reverse}{a logical stating whether the palette should be
    inverted (TRUE), or not (FALSE, default).}
  \item{x.min}{the minimal value from which to start the color scale}
  \item{x.max}{the maximal value from which to start the color scale}
  \item{na.col}{the color to be used for missing values (NAs)}
  \item{seed}{a seed for R's random number generated, used to fix the
    random permutation of colors in the palette used.}
  \item{col}{a vector of colors}
  \item{alpha}{a numeric value between 0 and 1 representing the alpha
  coefficient; 0: total transparency; 1: no transparency.}
}
\value{
  For \code{.genlab}, a character vector of size "n".
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## this opens the adegenet website
adegenetWeb()

}

.genlab("Locus-",11)

## transparent colors using "transp"
plot(rnorm(1000), rnorm(1000), col=transp("blue",.3), pch=20, cex=4)


## numeric values to color using num2col
plot(1:100, col=num2col(1:100), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=bluepal), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=flame), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=seasun), pch=20, cex=4)
plot(1:100, col=num2col(1:100, col.pal=seasun,rev=TRUE), pch=20, cex=4)

## factor as colors using fac2col
dat <- cbind(c(rnorm(50,8), rnorm(100), rnorm(150,3),
rnorm(50,10)),c(rnorm(50,1),rnorm(100),rnorm(150,3), rnorm(50,5)))
fac <- rep(letters[1:4], c(50,100,150,50))
plot(dat, col=fac2col(fac), pch=19, cex=4)
plot(dat, col=transp(fac2col(fac)), pch=19, cex=4)
plot(dat, col=transp(fac2col(fac,seed=2)), pch=19, cex=4)

}
\keyword{manip}