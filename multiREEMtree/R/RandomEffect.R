#' Return the random effect of a multiREEMtree object
#'
#' @param object A multiREEMtree object.
#'
#'
#' @return A list recording the random effects in a multiREEMtree object.
#'
#' @export

RandomEffect <- function(object){
   if(class(object)=="multiREEMtree"){
     return(object$RandomEffects)
   }else{
     stop("object should be a multiREEMtree.")
   }
}
