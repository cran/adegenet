\encoding{UTF-8}
\name{dapc}
\alias{dapc}
\alias{dapc.data.frame}
\alias{dapc.matrix}
\alias{dapc.genind}
\alias{print.dapc}
\alias{summary.dapc}
\alias{scatter.dapc}
\alias{assignplot}
\title{Discriminant Analysis of Principal Components (DAPC)}
\description{
  Important: this method is currently under review. Please email the
  author before using it.
  
  These functions implement the Discriminant Analysis of Principal Components
  (DAPC). See 'details' section for a succint description of the method. DAPC
  implementation calls upon \code{\link[ade4]{dudi.pca}} from the \code{ade4} package and
  \code{\link[MASS]{lda}} from the \code{MASS} package.

 \code{dapc} performs the DAPC on a \code{data.frame}, a \code{matrix}, or a
 \code{\linkS4class{genind}} object, and returns an object with class
 \code{dapc}. If data are stored in a \code{data.frame} or a \code{matrix},
 these have to be quantitative data (i.e., \code{numeric} or \code{integers}),
 as opposed to \code{characters} or \code{factors}.

 
 Other functions are:
  
  - \code{print.dapc}: prints the content of a \code{dapc} object.
  
  - \code{summary.dapc}: extracts useful information from a  \code{dapc} object.
  
  - \code{scatter.dapc}: produces scatterplots of principal components (or
    'discriminant functions'), with a screeplot of eigenvalues as inset.
  
  - \code{assignplot}: plot showing the probabilities of assignment of
    individuals to the different clusters.
}
\usage{
\method{dapc}{data.frame}(x, grp, n.pca=NULL, n.da=NULL, center=TRUE,
     scale=FALSE,var.contrib=FALSE, pca.select=c("nbEig","percVar"),
    perc.pca=NULL, \ldots)

\method{dapc}{matrix}(x, \ldots)

\method{dapc}{genind}(x, pop=NULL, n.pca=NULL, n.da=NULL, scale=FALSE,
     scale.method=c("sigma", "binom"), truenames=TRUE, all.contrib=FALSE,
     pca.select=c("nbEig","percVar"), perc.pca=NULL, \ldots)

\method{print}{dapc}(x, \dots)

\method{summary}{dapc}(object, \dots)

\method{scatter}{dapc}(x, xax=1, yax=2,
        col=rainbow(length(levels(x$grp))), posi="bottomleft", bg="grey",
        ratio=0.3, csub=1.2, \ldots)

assignplot(x, only.grp=NULL, subset=NULL, cex.lab=.75, pch=3)
}
\arguments{
\item{x}{\code{a data.frame}, \code{matrix}, or \code{\linkS4class{genind}}
  object. For the \code{data.frame} and \code{matrix} arguments, only
  quantitative variables should be provided.}
\item{grp,pop}{a \code{factor} indicating the group membership of individuals}
\item{n.pca}{an \code{integer} indicating the number of axes retained in the
  Principal Component Analysis (PCA) step. If \code{NULL}, interactive selection is triggered.}
\item{n.da}{an \code{integer} indicating the number of axes retained in the
  Discriminant Analysis step. If \code{NULL}, interactive selection is triggered.}
\item{center}{a \code{logical} indicating whether variables should be centred to
mean 0 (TRUE, default) or not (FALSE). Always TRUE for \linkS4class{genind} objects.}
\item{scale}{a \code{logical} indicating whether variables should be scaled
  (TRUE) or not (FALSE, default). Scaling consists in dividing variables by their
  (estimated) standard deviation to account for trivial differences in
  variances. Further scaling options are available for \linkS4class{genind}
  objects (see argument \code{scale.method}).}
\item{var.contrib,all.contrib}{a \code{logical} indicating whether the
  contribution of original variables (alleles, for \linkS4class{genind} objects)
  should be provided (TRUE) or not (FALSE, default). Such output can be useful,
  but can also create huge matrices when there is a lot of variables.}
\item{pca.select}{a \code{character} indicating the mode of selection of PCA
  axes, matching either "nbEig" or "percVar". For "nbEig", the user
  has to specify the number of axes retained (interactively, or via
  \code{n.pca}). For "percVar", the user has to specify the minimum amount of
  the total variance to be preserved by the retained axes, expressed as a
  percentage (interactively, or via \code{perc.pca}).  }
\item{perc.pca}{a \code{numeric} value between 0 and 100 indicating the
  minimal percentage of the total variance of the data to be expressed by the
  retained axes of PCA.}
\item{\ldots}{further arguments to be passed to other functions. For
  \code{dapc.matrix}, arguments are to match those of
  \code{dapc.data.frame}.}
\item{object}{a \code{dapc} object.}
\item{scale.method}{a \code{character} specifying the scaling method to be used
  for allele frequencies, which must match "sigma" (usual estimate of standard
  deviation) or "binom" (based on binomial distribution). See \code{\link{scaleGen}} for
  further details.}
\item{truenames}{a \code{logical} indicating whether true (i.e., user-specified)
  labels should be used in object outputs (TRUE, default) or not (FALSE).}
\item{xax,yax}{\code{integers} specifying which principal components of DAPC
  should be shown in x and y axes. }
\item{col}{a suitable color to be used for groups. The specified vector
should match the number of groups, not the number of individuals.}
\item{posi,bg,ratio,csub}{arguments used to customize the inset in scatterplots
  of DAPC results. See \code{\link[ade4]{add.scatter}} documentation in the
  ade4 package for
  more details.}
\item{only.grp}{a \code{character} vector indicating which groups should be
  displayed. Values should match values of \code{x$grp}. If \code{NULL}, all
  results are displayed}
\item{subset}{\code{integer} or \code{logical} vector indicating which
  individuals should be displayed. If \code{NULL}, all
  results are displayed}
\item{cex.lab}{a \code{numeric} indicating the size of labels.}
\item{pch}{a \code{numeric} indicating the type of point to be used to indicate
  the prior group of individuals (see \code{\link{points}} documentation for
  more details).}
}
\details{
  The Discriminant Analysis of Principal Components (DAPC) is designed
  to investigate the genetic structure of biological populations. This
  multivariate method consists in a two-steps procedure. First, genetic
  data are transformed (centred, possibly scaled) and submitted to a
  Principal Component Analysis (PCA). Second, principal components of
  PCA are submitted to a Linear Discriminant Analysis (LDA). A trivial
  matrix operation allows to express discriminant functions as linear
  combination of alleles, therefore allowing one to compute allele
  contributions. More details about the computation of DAPC are to be
  found in the indicated reference.

  DAPC does not infer genetic clusters ex nihilo; for this, see the
  \code{\link{find.clusters}} function.
}
\value{
  === dapc objects ===\cr
  The class \code{dapc} is a list with the following
  components:\cr
  \item{call}{the matched call.}
  \item{n.pca}{number of PCA axes retained}
  \item{n.da}{number of DA axes retained}
  \item{var}{proportion of variance conserved by PCA principal components}
  \item{eig}{a numeric vector of eigenvalues.}
  \item{grp}{a factor giving prior group assignment}
  \item{prior}{a numeric vector giving prior group probabilities}
  \item{assign}{a factor giving posterior group assignment}
  \item{tab}{matrix of retained principal components of PCA}
  \item{loadings}{principal axes of DAPC, giving coefficients of the linear
    combination of retained PCA axes.}
  \item{ind.coord}{principal components of DAPC, giving the coordinates of individuals onto
    principal axes of DAPC; also called the discriminant functions.}
  \item{grp.coord}{coordinates of the groups onto the principal axes of DAPC.}
  \item{posterior}{a data.frame giving posterior membership probabilities for
    all individuals and all clusters.}
  \item{var.contr}{(optional) a data.frame giving the contributions of original
    variables (alleles in the case of genetic data) to the principal components
    of DAPC.}


  === other outputs ===\cr
  Other functions have different outputs:\cr
  - \code{summary.dapc} returns a list with 6 components: \code{n.dim} (number
  of retained DAPC axes), \code{n.pop} (number of groups/populations),
  \code{assign.prop} (proportion of overall correct assignment),
  \code{assign.per.pop} (proportion of correct assignment per group),
  \code{prior.grp.size} (prior group sizes), and \code{post.grp.size} (posterior
  group sizes).

  - \code{scatter.dapc, assignplot} return the matched call.\cr
}
\references{
Jombart, T., Devillard, S. and Balloux, F.
Discriminant analysis of principal components: a new method for the analysis of
genetically structured populations. Submitted to \emph{BMC genetics}.
}
\seealso{
  - \code{\link{find.clusters}}: to identify clusters without prior.

  - \code{\link{dapcIllus}}: a set of simulated data illustrating the DAPC

  - \code{\link{eHGDP}}, \code{\link{H3N2}}: empirical datasets illustrating
  DAPC
}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
## data(dapcIllus), data(eHGDP), and data(H3N2) illustrate the dapc
## see ?dapcIllus, ?eHGDP, ?H3N2
##

example(dapcIllus)


\dontrun{
example(eHGDP)
example(H3N2)
}

}
\keyword{multivariate}