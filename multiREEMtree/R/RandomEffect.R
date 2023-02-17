#' Return the predicted random effects of a multiREEMtree or GmultiREEMtree object
#'
#'
#' @param object A multiREEMtree or GmultiREEMtree object.
#'
#'
#' @return A list recording the predicted random effects in a multiREEMtree or GmultiREEMtree object.
#'
#' @export

RandomEffect <- function(object){
   if(inherits(object, c("multiREEMtree", "GmultiREEMtree"))){
     return(object$RandomEffects)
   }else{
     stop("object shouldbe a multiREEMtree or GmultiREEMtree.")
   }
}
