\encoding{UTF-8}
\name{find.clusters}
\alias{find.clusters}
\alias{find.clusters.data.frame}
\alias{find.clusters.matrix}
\alias{find.clusters.genind}
\alias{.find.sub.clusters}
\title{find.cluster: cluster identification using successive K-means}
\description{
   These functions implement the clustering procedure used in Discriminant
  Analysis of Principal Components (DAPC, Jombart et al. submitted). This
  procedure consists in running successive K-means with an increasing number of
  clusters (\code{k}), after transforming data using a principal component
  analysis (PCA). For each model, a statistical measure of goodness of fit (by
  default, BIC) is computed, which allows to choose the optimal \code{k}. See
  \code{details} for a description of how to select the optimal \code{k}.

  Optionally, hierarchical clustering can be sought by providing a prior
  clustering of individuals (argument \code{clust}). In such case, clusters will
  be sought within each prior group.

  \code{.find.sub.clusters} is a hidden function called in some instances of
  \code{find.clusters}, and should not be called directly by the user.

  The K-means procedure used in \code{find.clusters} is \code{\link[stats]{kmeans}} function
  from the \code{stats} package. The PCA function is \code{\link[ade4]{dudi.pca}} from the
  \code{ade4} package.
}
\usage{
\method{find.clusters}{data.frame}(x, clust=NULL, n.pca=NULL,
              n.clust=NULL, stat=c("BIC","AIC", "WSS"),
              choose.n.clust=TRUE,criterion=c("min","diff", "conserv"),
              max.n.clust=round(nrow(x)/10), n.iter=1e3, n.start=10,
              center=TRUE, scale=TRUE, \ldots)

\method{find.clusters}{matrix}(x, \ldots)

\method{find.clusters}{genind}(x, clust=NULL, n.pca=NULL, n.clust=NULL,
              stat=c("BIC","AIC", "WSS"), choose.n.clust=TRUE,
              criterion=c("min","diff", "conserv"),
              max.n.clust=round(nrow(x@tab)/10), n.iter=1e3, n.start=10,
              scale=FALSE, scale.method=c("sigma", "binom"),
              truenames=TRUE, \ldots)

}
\arguments{
  \item{x}{\code{a data.frame}, \code{matrix}, or \code{\linkS4class{genind}}
    object. For the \code{data.frame} and \code{matrix} arguments, only
    quantitative variables should be provided.}
  \item{clust}{an optional \code{factor} indicating a prior group membership of
    individuals. If provided, sub-clusters will be sought within each prior
    group.}
  \item{n.pca}{an \code{integer} indicating the number of axes retained in the
    Principal Component Analysis (PCA) step. If \code{NULL}, interactive selection
    is triggered.}
  \item{n.clust}{ an optinal \code{integer} indicating the number of clusters to
    be sought. If provided, the function will only run K-means once, for this
    number of clusters. If left as \code{NULL}, several K-means are run for a
    range of k (number of clusters) values.}
  \item{stat}{ a \code{character} string matching 'BIC', 'AIC', or 'WSS', which
    indicates the statistic to be computed for each model (i.e., for each value of
    \code{k}). BIC: Bayesian Information Criterion. AIC: Aikaike's Information
    Criterion. WSS: within-groups sum of squares, that is, residual variance.}
  \item{choose.n.clust}{ a \code{logical} indicating whether the number of
    clusters should be chosen by the user (TRUE, default), or automatically,
    based on a given criterion (argument \code{criterion}). IT IS HIGHLY
    RECOMMENDED to choose the number of clusters interactively, as automatic
    procedures have not been fully evaluated.}
  \item{criterion}{ a \code{character} string matching "min", "diff", or
    "conserv", indicating the criterion for automatic selection of the optimal
    number of clusters. Honestly, you should go for interactive
    selection of the number of clusters. Do as you wish. No warranty. If
    you still want to give it a try, see \code{details}.}
  \item{max.n.clust}{ an \code{integer} indicating the maximum number of
    clusters to be tried. Values of 'k' will be picked up between 1 and \code{max.n.clust}}
  \item{n.iter}{ an \code{integer} indicating the number of iterations to be used
    in each run of K-means algorithm. Corresponds to \code{iter.max} of
    \code{kmeans} function.}
  \item{n.start}{ an \code{integer} indicating the number of randomly
    chosen starting centroids to be used in each run of the K-means
    algorithm. Using more starting points ensures convergence of the
    algorithm. Corresponds to \code{nstart} of \code{kmeans} function.}
  \item{center}{a \code{logical} indicating whether variables should be centred to
    mean 0 (TRUE, default) or not (FALSE). Always TRUE for \linkS4class{genind}
    objects.}
  \item{scale}{a \code{logical} indicating whether variables should be
    scaled (TRUE) or not (FALSE, default). Scaling consists in dividing
    variables by their (estimated) standard deviation to account for
    trivial differences in variances. In allele frequencies, it comes with
    the risk of giving uninformative alleles more importance while
    downweighting informative alleles. Further scaling options are
    available for \linkS4class{genind} objects (see argument
    \code{scale.method}).}
  \item{scale.method}{a \code{character} specifying the scaling method to be used
    for allele frequencies, which must match "sigma" (usual estimate of standard
    deviation) or "binom" (based on binomial distribution). See
    \code{\link{scaleGen}} for further details.}
  \item{truenames}{a \code{logical} indicating whether true (i.e., user-specified)
    labels should be used in object outputs (TRUE, default) or not
    (FALSE), in which case generic labels are used.}
  \item{\ldots}{further arguments to be passed to other functions. For
    \code{find.clusters.matrix}, arguments are to match those of the
    \code{data.frame} method.}
}
\details{
  === ON THE SELECTION OF K ===\cr
  (where K is the 'optimal' number of clusters)

  So far, the analysis of data simulated under various population genetics
  models (see reference) suggested an ad hoc rule for the selection of the
  optimal number of clusters. First important result is that BIC seems for
  efficient than AIC and WSS to select the appropriate number of clusters (see
  example). The rule of thumb consists in increasing K until it no longer leads
  to an appreciable improvement of fit (i.e., to a decrease of BIC).  In the
  most simple models (island models), BIC decreases until it reaches the optimal
  K, and then increases. In these cases, our rule amounts to choosing the lowest
  K. In other models such as stepping stones, the decrease of BIC often
  continues after the optimal K, but is much less steep.

  
  An alternative approach, that we do not recommend for now, is automatic
  selection based on a fixed criterion. For this, set \code{choose.n.clust} to
  FALSE and specify the \code{criterion} you want to use, from the following
  values:

  - "min": the model with the minimum summary statistics (as specified by
    \code{stat} argument, BIC by default) is retained. Is likely to work for
    simple island model, using BIC.

  - "diff": model selection based on successive improvement of the test
  statistic. This procedure attempts to increase K until the model improvement
  (difference in successive BIC, AIC, or WSS) is no longer important. May be
  more appropriate to models relating to stepping stones.

  - "conserv": another criterion meant to be conservative, in that it seeks a
  good fit with a minimum number of clusters. Unlike "diff", it does not rely on
  differences between successive statistics, but rather on absolute fit. It
  selects the model with the smallest K so that the overall fit is above a given
  threshold.
}
\value{
  The class \code{find.clusters} is a list with the following
  components:\cr
  \item{Kstat}{a \code{numeric} vector giving the values of the summary
  statistics for the different values of K. Is NULL if \code{n.clust} was specified.}
  \item{stat}{a \code{numeric} value giving the value of the summary statistics
  for the retained model}
  \item{grp}{a \code{factor} giving group membership for each individual.}
  \item{size}{an \code{integer} vector giving the size of the different clusters.}
}
\references{
  Jombart T, Devillard S and Balloux F  (2010) Discriminant analysis of
  principal components: a new method for the analysis of genetically
  structured populations. BMC Genetics
  11:94. doi:10.1186/1471-2156-11-94
}
\seealso{
  - \code{\link{dapc}}: implements the DAPC.

  - \code{\link{dapcIllus}}: dataset illustrating the DAPC and \code{find.clusters}.

  - \code{\link{eHGDP}}: dataset illustrating the DAPC and \code{find.clusters}.

  - \code{\link[stats]{kmeans}}: implementation of K-means in the stat package.

  - \code{\link[ade4]{dudi.pca}}: implementation of PCA in the ade4 package.

}
\author{ Thibaut Jombart \email{t.jombart@imperial.ac.uk} }
\examples{
\dontrun{
## THIS ONE TAKES A FEW MINUTES TO RUN ## 
data(eHGDP)

## here, n.clust is specified, so that only on K value is used
grp <- find.clusters(eHGDP, max.n=30, n.pca=200, scale=FALSE, n.clust=4) # takes about 2 minutes
names(grp)
grp$Kstat
grp$stat


## to try different values of k (interactive)
grp <- find.clusters(hgdp, max.n=50, n.pca=200, scale=FALSE)

## and then, to plot BIC values:
plot(grp$Kstat, type="b", col="blue")
}


## ANOTHER SIMPLE EXAMPLE ## 
data(sim2pop) # this actually contains 2 pop

## DETECTION WITH BIC (clear result)
foo.BIC <- find.clusters(sim2pop, n.pca=100, choose=FALSE)
plot(foo.BIC$Kstat, type="o", xlab="number of clusters (K)", ylab="BIC",
col="blue", main="Detection based on BIC")
points(2, foo.BIC$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")


## DETECTION WITH AIC (less clear-cut)
foo.AIC <- find.clusters(sim2pop, n.pca=100, choose=FALSE, stat="AIC")
plot(foo.AIC$Kstat, type="o", xlab="number of clusters (K)", ylab="AIC", col="purple", main="Detection based on AIC")
points(2, foo.AIC$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")


## DETECTION WITH WSS (less clear-cut)
foo.WSS <- find.clusters(sim2pop, n.pca=100, choose=FALSE, stat="WSS")
plot(foo.WSS$Kstat, type="o", xlab="number of clusters (K)", ylab="WSS
(residual variance)", col="red", main="Detection based on WSS")
points(2, foo.WSS$Kstat[2], pch="x", cex=3)
mtext(3, tex="'X' indicates the actual number of clusters")

}
\keyword{multivariate}
