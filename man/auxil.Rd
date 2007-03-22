\name{Auxiliary functions}
\alias{truenames}
\alias{.truenames.genind}
\alias{.truenames.genpop}
\alias{seploc}
\alias{.seploc.genind}
\alias{.seploc.genpop}
\title{ Utilities functions for adegenet objects}
\description{ Theses functions are to be used with \code{genind} and
  \code{genpop} objects.\cr
  \code{truenames} returns elements of the object with using true names
  (as opposed to generic labels) for individuals, markers, alleles, and population.\cr
  The function \code{seploc} splits the table (\code{x\$tab}) by marker,
  allowing separate analysis of markers.
}
\usage{
truenames(x)
seploc(x,truenames=FALSE)
}
\arguments{
  \item{x}{ a \code{genind} or \code{genpop} object.}
  \item{truenames}{a logical indicating whether generic labels
    (FALSE,default) or true names should be used (TRUE).}
}
\value{The function \code{truenames} returns a matrix similar to x\$tab
but with true labels. If x\$pop exists, it returns a list with this
matrix (\$tab) and a population vector with true names (\$pop).\cr
   
   The function \code{seploc} applied to \code{genind} or \code{genpop}
objects returns a list of matrices, one per marker.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
data(microbov)
# restore true names
truenames(microbov)$tab[1:5,1:5]

# isolate each marker
obj <- seploc(microbov,truenames=TRUE)
names(obj)

# make a new object with INRA5
head(obj$INRA5)
inra5.gind <- as.genind(obj$INRA5)
inra5.gind

# perform tests only on this marker
if(require(genetics)){
hw.test <- HWE.test.genind(inra5.gind,pop=microbov$pop,res.type="matrix",permut=TRUE) 
hw.test
}

if(require(hierfstat)){
g.test <- gstat.randtest(inra5.gind,pop=microbov$pop,nsim=99)
g.test
}
}
\keyword{manip}