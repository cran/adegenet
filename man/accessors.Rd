\encoding{UTF-8}
\docType{methods}
\name{Accessors}
\alias{$,genind-method}
\alias{$,genpop-method}
\alias{$<-,genind-method}
\alias{$<-,genpop-method}
\alias{[,genind-method}
\alias{[,genpop-method}
\alias{nLoc}
\alias{nLoc,genind-method}
\alias{nLoc,genpop-method}
\alias{pop}
\alias{pop<-}
\alias{pop,genind-method}
\alias{pop<-,genind-method}
\alias{locNames}
\alias{locNames,genind-method}
\alias{locNames,genpop-method}
\title{ Accessors for adegenet objects}
\description{
  An accessor is a function that allows to interact with slots of an
  object in a convenient way. Several accessors are available for \linkS4class{genind} or
  \linkS4class{genpop} objects. The operator "\$" and "\$<-" are used to
  access the slots, being equivalent to "@" and "@<-".\cr

  The operator "[" can be used to access components of the matrix slot
  "@tab", returning a \linkS4class{genind} or \linkS4class{genpop}
  object. This syntax is the same as for a matrix; for instance:\cr
  - "obj[,]" returns "obj" \cr
  - "obj[1:10,]" returns an object with only the first 10 genotypes (if
  "obj" is a \linkS4class{genind}) or the first 10 populations (if
  "obj" is a \linkS4class{genpop}) of "obj" \cr
  - "obj[1:10, 5:10]" returns an object keeping the first 10 entities and
  the alleles 5 to 10.\cr
  - "obj[loc=c("L1","L3")]" returns an object keeping only the loci
  specified in the \code{loc} argument (using generic names, not true
  names; in this example, only the first and the third locus would be
  retained)\cr
  - "obj[1:3, drop=TRUE]" returns the first 3 genotypes/populations of
  "obj", but retaining only alleles that are present in this subset (as
  opposed to keeping all alleles of "obj", which is the default
  behavior).\cr

  The argument \code{treatOther} handles the treatment of objects in the
  \code{@other} slot (see details). The argument \code{drop} can be set
  to TRUE to drop alleles that are no longer represented in the subset.
}
\section{Methods}{
  \describe{
    \item{nLoc}{returns the number of loci of the object}
    \item{pop}{returns the population factor of the object, using true
      (as opposed to generic) levels.}
    \item{pop<-}{replacement method for the \code{@pop} slot of an
      object. The content of \code{@pop} and \code{@pop.names} is updated
      automatically.}
    \item{locNames}{returns the true names of markers and/or alleles.}
  }
}
\usage{
nLoc(x, \dots)
pop(x)
locNames(x, \dots)
\S4method{locNames}{genind}(x, withAlleles=FALSE, \dots)
\S4method{locNames}{genpop}(x, withAlleles=FALSE, \dots)
}
\arguments{
  \item{x}{a \linkS4class{genind} or a \linkS4class{genpop} object.}
  \item{withAlleles}{a logical indicating whether the result should be
    of the form [locus name].[allele name], instead of [locus name].}
  \item{\dots}{further arguments to be passed to other methods
  (currently not used).}
}
\value{
  A \linkS4class{genind} or \linkS4class{genpop} object.
}
\author{Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\details{
  The "[" operator can treat elements in the \code{@other} slot as
  well. For instance, if \code{obj@other$xy} contains spatial
  coordinates, the \code{obj[1:3,]@other$xy} will contain the spatial
  coordinates of the genotypes (or population) 1,2 and 3. This is
  handled through the argument \code{treatOther}, a logical defaulting
  to TRUE. If set to FALSE, the \code{@other} component is not
  returned.\cr
  
  Note that only matrix-like, vector-like and lists can be proceeded in
  \code{@other}. Other kind of objects will issue a warning an be
  returned as they are.\cr

  The \code{drop} argument can be set to TRUE to retain only alleles
  that are present in the subset. To achieve better control of
  polymorphism of the data, see \code{\link{isPoly}}.
}
\examples{
data(nancycats)
nancycats
nancycats$pop

# let's isolate populations 4 and 8
temp <- nancycats@pop=="P04" | nancycats@pop=="P08"
obj <- nancycats[temp,]
obj

truenames(obj)$pop

# let's isolate two markers, fca23 and fca90
nancycats$loc.names

# they correspond to L2 and L7
temp <- nancycats$loc.fac=="L2" | nancycats$loc.fac=="L7"
obj <- nancycats[,temp]
obj

obj$loc.fac 
obj$loc.names

# or more simply
nancycats[loc=c("L2","L7")]
obj$loc.fac 
obj$loc.names

# using 'drop':
truenames(nancycats[1:2])$tab
truenames(nancycats[1:2, drop=TRUE])$tab

# illustrate how 'other' slot is handled
colonies <- genind2genpop(nancycats)
colonies@other$aChar <- "This will not be proceeded"
colonies123 <- colonies[1:3]
colonies
colonies@other$xy

# illustrate pop
obj <- nancycats[sample(1:100,10)]
obj$pop
obj$pop.names
pop(obj)
pop(obj) <- rep(c('b','a'), each=5)
obj$pop
obj$pop.names
pop(obj)

# illustrate locNames
locNames(obj)
locNames(obj, withAlleles=TRUE)

}
\keyword{manip}