#' Plot the multivariate tree of a GmultiREEMtree object
#'
#' Plot the multivariate tree of a GmultiREEMtree object.
#'
#' @param x A GmultiREEMtree object.
#' @param ... Identical to the arguments in \link[mvpart]{text.rpart}
#'
#' @return A plot of a multivariate RE-EM tree.
#'
#' @import mvpart
#' @import MASS
#' @importFrom graphics text
#' @rdname plot.GmultiREEMtree
#' @export

plot.GmultiREEMtree<- function(x, ...){
   if(inherits(x, "GmultiREEMtree")){

     plot_mvpart <- utils::getS3method('plot', 'rpart', envir = asNamespace("mvpart"))

     plot_mvpart(x$Tree)

     text_mvpart <- utils::getS3method('text', 'rpart', envir = asNamespace("mvpart"))

     text_mvpart(x$Tree, ...)
   }
}
