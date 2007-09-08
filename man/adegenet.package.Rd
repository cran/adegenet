\encoding{UTF-8}
\name{adegenet-package}
\alias{adegenet-package}
\alias{adegenet}
\docType{package}
\title{Genetic data handling for multivariate analysis}

\description{This package is devoted to manipulate data obtained from
  molecular markers. The newly defined classes of object facilitate their
  analysis within the multivariate framework of the \code{ade4}
  package. However, this package also provides interfaces with other
  packages, making Hardy-Weinberg equilibrium test, F statistics,
  Goudet's G test for population structure, or linkage disequilibrium
  measure available directly or using simple conversion functions.\cr
  Since the second version (1.0-1), \code{adegenet} also includes
  spatial analysis function like Monmonier algorithm (see
  \code{\link{monmonier}}) to define genetic boundaries among
  individuals or populations.

  The basic class of object is \code{\link{genind}}, and contains
  genotypes (genind stands for genotypes-individuals).\cr
  It can be obtained by converting files from GENETIX, Fstat
  and Genepop using \code{\link{import2genind}}.\cr
  The second class is \code{\link{genpop}}: such object contains alleles counts per populations and
  loci. It can be obtained from any \code{\link{genind}} object using
  \code{\link{genind2genpop}}.\cr
  It is also possible to obtain a table of allelic frequencies using
  \code{\link{makefreq}} on a \code{\link{genpop}} object.
  In all cases, missing data can be treated using different options.\cr

  The package proposes useful functions for \code{genind} objects:\cr
  - \code{\link{HWE.test.genind}} to test for Hardy-Weinberg equilibrium on every locus
  x population combinaisons (based on \code{HWE.test}, package
  \code{genetics}).\cr
  - \code{\link{gstat.randtest}} is Monte Carlo test (class
  \code{\link[pkg:ade4]{randtest}}) of Goudet's G statistic measuring
  population structure (based on \code{\link[pkg:hierfstat]{g.stat.glob}}, package
  \code{hierfstat}).\cr
  - \code{\link{genind2genotype}} assures conversion into
  \code{\link[pkg:genetics]{genotype}} objects used in \code{genetics}
  and \code{LDheatmap} packages.\cr
  - \code{\link{genind2hierfstat}} assures conversion into
  the format (particular data.frame) used in \code{hierfstat} package.\cr

  Moreover, several genetic distances between populations can be computed
  using \code{\link{dist.genpop}}.
  Genetic distances between individuals are not yet implemented.

  Lastly, the Monmonier algorithm is implemented to seek genetic
  boundaries. The main function is \code{\link{monmonier}}, but
  better (optimized) boundaries can be obtained using
  \code{\link{optimize.monmonier}}. Object of the class \code{monmonier}
  can be plotted and printed using the corresponding methods.   
  
}
\details{
  \tabular{ll}{
    Package: \tab adegenet\cr
    Type: \tab Package\cr
    Version: \tab 1.0-2\cr
    Date: \tab 2007-09-08\cr
    License: \tab GPL version 2 or newer
  }  

  These are the essential functions provided by the package:\cr

  \code{\link{genind}}: adegenet class for individual genotypes\cr
  
  \code{\link{genpop}}: adegenet class for allele counts in
  populations\cr
  
  \code{\link{import2genind}}: Conversion function for adegenet (from
  GENETIX, Fstat, Genepop)\cr
  
  \code{\link{makefreq}}: Function to generate allelic frequencies\cr

  \code{\link{HWE.test.genind}}: Hardy-Weinberg Equilibrium test for
  multilocus data\cr

  \code{\link{gstat.randtest}}: Monte Carlo test of Goudet's G statistic
  for multilocus data\cr

  \code{\link{monmonier}}: Monmonier algorithm to find genetic
  boundaries among genetic entities on a connection network\cr

}
\author{Thibaut Jombart <jombart@biomserv.univ-lyon1.fr>}
\references{
  See adegenet website: \url{http://pbil.univ-lyon1.fr/software/adegenet/}
}
\keyword{manip}
\keyword{multivariate}
\seealso{\code{ade4} package for multivariate analysis}