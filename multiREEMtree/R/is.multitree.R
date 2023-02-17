#' Check if an R object is a multitree object
#'
#' Check if an R object is a multitree object.
#'
#' @param object An object to be checked.
#'
#' @return A logical value
#'
#' @import mvpart
#' 
#' @rdname is.multitree
#' @export

is.multitree <- function(object){
  return(inherits(object, "multitree"))
}
