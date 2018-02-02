#' snapclust's Akaike Information Criterion (AIC)
#'
#' This function computes Akaike Information Criterion (AIC) for
#' \code{snapclust} results.
#'
#' @author Thibaut Jombart \email{thibautjombart@@gmail.com}
#'
#' @export
#'
#' @param object An object returned by the function \code{\link{snapclust}}.
#'
#' @param ... Further arguments for compatibility with the \code{AIC} generic
#'     (currently not used).
#'
#' @references Beugin M-P, Gayet T, Pontier D, Devillard S, Jombart T. A fast
#'     likelihood solution to the genetic clustering problem. Methods Ecol
#'     Evol. 2018;00:1–11. \url{https://doi.org/10.1111/2041-210X.12968}
#'
#' @seealso
#' \itemize{
#'  \item \code{\link{snapclust}}: to identify clusters
#'
#'  \item \code{\link{snapclust.choose.k}}: to find the number of clusters
#'
#'  \item \code{\link{AICc.snapclust}}: AICc computation
#'
#'  \item \code{\link{BIC.snapclust}}: BIC computation
#'
#'  \item \code{\link{KIC.snapclust}}: KIC computation
#' }
#'
#'
AIC.snapclust <- function(object, ...) {

    ## The number of parameters is defined as:
    ## (number of independent allele frequencies) x (nb clusters).

  -2 * object$ll + 2 * object$n.param
}
