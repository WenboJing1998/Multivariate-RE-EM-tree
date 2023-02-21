#' Return the fitted values of a multitree object
#'
#' Return the fitted values of a multitree object.
#'
#' @param object A multitree object.
#' @param ... Won't be used.
#'
#' @return A matrix of fitted values.
#'
#' @import mvpart
#' 
#' @rdname fitted.multitree
#' @export


fitted.multitree <- function(object, ...){

    if(inherits(object, "multitree")){

      return(object$fitted)

    }else{
      NextMethod("fitted")
    }
}

