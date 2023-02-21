#' Check if an R object is a GmultiREEMtree object
#'
#' Check if an R object is a GmultiREEMtree object.
#'
#' @param object An object to be checked.
#'
#' @return A logical value
#'
#' @import mvpart
#' 
#' @rdname is.GmultiREEMtree
#' @export

is.GmultiREEMtree <- function(object){
  return(inherits(object, "GmultiREEMtree"))
}
