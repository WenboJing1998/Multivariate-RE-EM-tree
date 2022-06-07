#' Summarize a multiREEMtree object
#'
#' Summarize a multiREEMtree object.
#'
#'
#' @param object a multiREEMtree object.
#' @param ... Won't be used.
#'
#'
#' @rdname summary.multiREEMtree
#' @export

summary.multiREEMtree <- function(object,...){

   if(class(object)=="multiREEMtree"){
     print("Multivariate Tree:")
     summary(object$Tree)
     print("Mixed Effect Model:")
     summary(object$EffectModel)

    } else{
      NextMethod("summary")
    }
  }

