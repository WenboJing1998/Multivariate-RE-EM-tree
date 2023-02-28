#' Return the fitted values of a multiREEMtree object
#'
#' Return the fitted values of a multiREEMtree object.
#'
#' @param object A multiREEMtree object.
#' @param level  An optional integer giving the level of grouping to be used in
#' extracting the fitted values from object. Defaults to be the innermost level (level=-1).
#' If the lme.method of the object is "MCMCglmm", fitted levels can only be 0 (for fixed) or 1
#' (for fixed + random).
#' @param ... Won't be used.
#'
#' @return A matrix of fitted values.
#'
#' @import mvpart
#' @import MASS
#' @importFrom stats fitted
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
#' result <- multiREEMtree(formula=f, data=data[1:90, ], random=r)
#' fit <- fitted(result)
#' @rdname fitted.multiREEMtree
#' @export


fitted.multiREEMtree <- function(object, level=-1, ...){

    if(inherits(object, "multiREEMtree")){

      maxlevel <-  length(object$RandomEffects)

      if(level > maxlevel) stop(paste("The max level for this object is", maxlevel))

      if(level==-1) level <- maxlevel

      return(object$fitted[[level+1]])

    }else{
      NextMethod("fitted")
    }
}

