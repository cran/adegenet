\encoding{UTF-8}
\name{genpop}
\alias{genpop}
\alias{is.genpop}
\alias{as.genpop}
\alias{print.genpop}
\alias{summary.genpop}
\title{adegenet class for allele counts in populations}
\description{The objects of class \code{genpop} contain alleles counts
  for several loci.\cr
  It consists in a list with several components (see value section).\cr
  Such object is obtained using \code{genind2genpop} which converts
  individuals genotypes of known population into a \code{genpop} object.
  Note that the function \code{summary} of a \code{genpop} object
  returns a list of components.
}
\usage{
is.genpop(x)
as.genpop(tab = NULL, prevcall = NULL)
\method{print}{genpop}(x, \dots)
\method{summary}{genpop}(object, \dots)
}
\arguments{
  \item{x}{an object of class \code{genpop}.}
  \item{tab}{a populations x alleles matrix of allele counts.}
  \item{prevcall}{call of an object, for internal use.}
  \item{...}{other -unused- arguments}
  \item{object}{an object of class \code{genpop}.}
}
\value{
  \item{tab}{matrix of alleles counts for each combinaison of population
    -in rows- and alleles -in columns-. Rows and columns are given generic names.}
  \item{pop.names}{character vector containing the real names of the populations}
  \item{loc.names}{character vector containing the real names of the loci}
  \item{loc.nall}{integer vector giving the number of alleles per locus}
  \item{loc.fac}{locus factor for the columns of \code{tab}}
  \item{all.names}{list having one component per locus, each containing a character vector of alleles names}
  \item{call}{the matched call}
  \item{npop}{(summary) number of populations.}
  \item{loc.nall}{(summary) number of alleles per locus.}
  \item{pop.nall}{(summary) number of alleles per population.}
  \item{NA.perc}{(summary) percentage of - appearing - missing data.}
}
\references{}
\seealso{\code{\link{makefreq}}, \code{\link{genind}}, \code{\link{import2genind}}, \code{\link{genetix2genind}}, 
\code{\link{genepop2genind}}, \code{\link{fstat2genind}}
}
\author{ Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
\examples{
obj1 <- import2genind(system.file("files/nancycats.gen",
package="adegenet"))
is.genpop(obj1)
summary(obj1)
obj1


obj2 <- genind2genpop(obj1)
is.genpop(obj2)
obj2

if(require(ade4)){
data(microsatt)
# use as.genpop to convert convenient count tab to genpop
obj3 <- as.genpop(microsatt$tab)
obj3

all(obj3$tab==microsatt$tab)
all(obj3$pop.names==rownames(microsatt$tab))
# it worked

# perform a correspondance analysis
obj4 <- genind2genpop(obj1,missing="replace")
ca1 <- dudi.coa(as.data.frame(obj4$tab),scannf=FALSE)
s.label(ca1$li,sub="Correspondance Analysis",csub=2)
add.scatter.eig(ca1$eig,2,xax=1,yax=2,posi="top")
}
}
\keyword{manip}
\keyword{multivariate}

