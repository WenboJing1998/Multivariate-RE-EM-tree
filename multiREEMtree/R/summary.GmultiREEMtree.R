#' Summarize a GmultiREEMtree object
#'
#' Summarize a GmultiREEMtree object.
#'
#'
#' @param object a GmultiREEMtree object.
#' @param ... Won't be used.
#'
#'
#' @rdname summary.GmultiREEMtree
#' @export

summary.GmultiREEMtree <- function(object,...){

   if(inherits(object, "GmultiREEMtree")){
     print("Multivariate Tree:")
     summary(object$Tree)
     print("Mixed Effect Model:")
     summary(object$EffectModel)

    } else{
      NextMethod("summary")
    }
  }

