#' Plot the multivariate tree of a multitree object
#'
#' Plot the multivariate tree of a multitree object.
#'
#' @param x A multitree object.
#' @param ... Identical to the arguments in \link[mvpart]{text.rpart}
#'
#' @return A plot of a multivariate RE-EM tree.
#'
#' @import mvpart
#' @import MASS
#' @importFrom graphics text
#' @examples
#' require(multiREEMtree)
#' # simulate multivariate grouped data
#' X1 <- rnorm(100, 0, 1)
#' X2 <- rnorm(100, 0, 1)
#' X3 <- rnorm(100, 0, 1)
#' group <- rep(1:10, 10)
#' simb <- MASS::mvrnorm(10, rep(0, 2), diag(2))
#' sime <- MASS::mvrnorm(100, rep(0, 2), diag(2))
#' Y1 <- 10 * (X1>0) + 20 * (X1<0 & X2 >0) + simb[group, 1] + sime[ ,1]
#' Y2 <- 30 * (X1>0) + 50 * (X1<0 & X2 >0) + simb[group, 2] + sime[ ,2]
#' data <- data.frame(cbind(X1, X2, X3, Y1, Y2, group))
#' f <- Y1 + Y2 ~ X1 + X2 + X3
#' result <- multitree(formula=f, data=data)
#' plot(result)
#' @rdname plot.multitree
#' @export

plot.multitree <- function(x, ...){
   if(inherits(x, "multitree")){

     plot_mvpart <- utils::getS3method('plot', 'rpart', envir = asNamespace("mvpart"))

     plot_mvpart(x$Tree)

     text_mvpart <- utils::getS3method('text', 'rpart', envir = asNamespace("mvpart"))

     text_mvpart(x$Tree, ...)
   }
}
