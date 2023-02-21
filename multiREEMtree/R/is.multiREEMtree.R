#' Check if an R object is a multiREEMtree object
#'
#' Check if an R object is a multiREEMtree object.
#'
#' @param object An object to be checked.
#'
#' @return A logical value
#'
#' @import mvpart
#' @import MASS
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
#' Y2 <- 30 * (X1>0) + 50 * (X1<0 & X2 >0)+ simb[group, 2] + sime[ ,2]
#' data <- data.frame(cbind(X1, X2, X3, Y1, Y2, group))
#' f <- Y1 + Y2 ~ X1 + X2 + X3
#' r <-  ~ 1 | group
#' result <- multiREEMtree(formula=f, data=data, random=r)
#' is.multiREEMtree(result)
#' @rdname is.multiREEMtree
#' @export

is.multiREEMtree <- function(object){
  return(inherits(object, "multiREEMtree"))
}
