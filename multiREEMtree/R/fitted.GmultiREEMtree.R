#' Return the fitted values of a GmultiREEMtree object
#'
#' Return the fitted values of a GmultiREEMtree object.
#'
#' @param object A GmultiREEMtree object.
#' @param level  An optional integer giving the level of grouping to be used in
#' extracting the fitted values from object. Defaults to be the innermost level (level=-1).
#' @param ... Won't be used.
#'
#' @return A matrix of fitted values.
#'
#' @import mvpart
#' @importFrom stats fitted
#' 
#' @rdname fitted.GmultiREEMtree
#' @export


fitted.GmultiREEMtree <- function(object, level=-1, ...){

    if(inherits(object, "GmultiREEMtree")){

      maxlevel <-  length(object$RandomEffects)

      if(level > maxlevel) stop(paste("The max level for this object is", maxlevel))

      if(level==-1) level <- maxlevel

      return(object$fitted[[level+1]])

    }else{
      NextMethod("fitted")
    }
}

