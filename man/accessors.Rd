\encoding{UTF-8}
\docType{methods}
\name{Accessors}
\alias{$,genind-method}
\alias{$,genpop-method}
\alias{$<-,genind-method}
\alias{$<-,genpop-method}
\alias{[,genind-method}
\alias{[,genpop-method}
\title{ Accessors for adegenet objects}
\description{
  Several accessors for \linkS4class{genind} or
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
  - "obj[nloc=c("L1","L3")]" returns an object keeping only the loci
  specified in the \code{loc} argument (using generic names, not true
  names; in this example, only the first and the third locus would be
  retained)\cr
  
  The argument \code{treatOther} handles the treatment of objects in the
  \code{@other} slot (see details).
}
\usage{
}
\value{
  A \linkS4class{genind} or \linkS4class{genpop} object.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
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

# illustrate how 'other' slot is handled
colonies <- genind2genpop(nancycats)
colonies@other$aChar <- "This will not be proceeded"
colonies123 <- colonies[1:3]
colonies
colonies@other$xy
}
\keyword{manip}