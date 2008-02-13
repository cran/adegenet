\encoding{UTF-8}
\docType{methods}
\name{Accessors}
\alias{$,pdata-method}
\alias{$,phylo4-method}
\alias{$<-,pdata-method}
\alias{$<-,phylo4-method}
\alias{[,genind-method}
\alias{[,genpop-method}
\title{ Accessors for adegenet objects}
\description{
  Several accessors for \linkS4class{genind} or
  \linkS4class{genpop} objects. The operator "\$" and "\$<-" are used to
  access the slots, being equivalent to "@" and "@<-".

  The operator "[" can be used to access components of the matrix slot
  "@tab", returning a \linkS4class{genind} or \linkS4class{genpop}
  object. This syntax is the same as for a matrix; for instance:\cr
  - "obj[,]" returns "obj" \cr
  - "obj[1:10,]" returns an object with only the first 10 genotypes (if
  "obj" is a \linkS4class{genind}) or the first 10 populations (if
  "obj" is a \linkS4class{genpop}) of "obj" \cr
  - "obj[1:10, 5:10]" returns an object keeping the first 10 entities and
  the alleles 5 to 10. 
  
}
\usage{
}
\value{
  A \linkS4class{genind} or \linkS4class{genpop} object.
}
\author{Thibaut Jombart \email{jombart@biomserv.univ-lyon1.fr} }
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
}
\keyword{manip}