#' Return the residuals of a multitree object
#'
#' @param object A multitree object.
#' @param ... Won't be used.
#'
#' @return A matrix of residuals.
#'
#' @import mvpart
#' 
#' @rdname residuals.multitree
#' @export


residuals.multitree <- function(object, ...){

    if(inherits(object, "multitree")){

        return(object$residuals)

    }else{
      NextMethod("residuals")
    }
  }
