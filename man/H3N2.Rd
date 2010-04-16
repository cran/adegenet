\encoding{UTF-8}
\name{H3N2}
\alias{H3N2}
\docType{data}
\title{Seasonal influenza (H3N2) HA segment data}
\description{
This dataset consists of 1903 strains of seasonal influenza (H3N2)
distributed worldwide, and typed at 125 SNPs located in the
hemagglutinin (HA) segment. These data were gathered from DNA sequences
available from Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).
}
\usage{data(H3N2)}
\format{
    \code{H3N2} is a genind object with several data frame as supplementary
    components (\code{H3N2@other) slort}, which contains the following items:
    \describe{
        \item{x}{a \code{data.frame} containing miscellanous annotations
	  of the sequences.}
        \item{xy}{a matrix with two columns indicating the geographic
	  coordinates of the strains, as longitudes and latitudes.}
        \item{epid}{a character vector indicating the epidemic of the strains.}
    }
}
\source{
  This dataset was prepared by Thibaut Jombart
  (t.jombart@imperia.ac.uk), from annotated sequences available on
  Genbank (http://www.ncbi.nlm.nih.gov/Genbank/).
}
\references{
Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{Genetics}.
}
\examples{
## LOAD DATA
data(H3N2)
H3N2

## set population to yearly epidemics
pop(H3N2) <- factor(H3N2$other$epid)



## PERFORM DAPC - USE POPULATIONS AS CLUSTERS
## to reproduce exactly analyses from the paper, use "n.pca=1000"
dapc1 <- dapc(H3N2, all.contrib=TRUE, scale=FALSE, n.pca=150, n.da=5)
dapc1

## (see ?dapc for details about the output)


## SCREEPLOT OF EIGENVALUES
barplot(dapc1$eig, main="H3N2 - DAPC eigenvalues")


## SCATTERPLOT (axes 1-2)
scatter(dapc1, ratio=.2)

}
\keyword{datasets}
