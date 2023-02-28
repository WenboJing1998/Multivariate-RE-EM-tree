#' Summarize a multitree object
#'
#' Summarize a multitree object.
#'
#'
#' @param object a multiree object.
#' @param ... Won't be used.
#'
#'
#' @rdname summary.multitree
#' @export

summary.multitree <- function(object,...){

   if(inherits(object, "multitree")){
     print("Multivariate Tree:")
     summary(object$Tree)
    } else{
      NextMethod("summary")
    }
  }

