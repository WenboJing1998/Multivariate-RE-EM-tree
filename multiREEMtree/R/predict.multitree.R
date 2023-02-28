#' Make predictions based on a multitree object on a new dataset
#'
#' Make predictions based on a multitree object on a new dataset.
#'
#' @param object A multiree object.
#' @param newdata A data frame in which to look for variables with which to predict.
#' @param ... Won't be used.
#' @return A matrix of predicted values.
#'
#' @import mvpart
#' @import MASS
#' @import expm
#' @import utils
#' @importFrom stats predict sd cov
#' 
#' @rdname predict.multitree
#' @export

predict.multitree <- function(object, newdata, ...){
  if(inherits(object, "multitree")){
    
    VariableNumber <- length(object$ResponseVariables)
    
    TreePrediction <- matrix(NA, nrow = dim(newdata)[1], ncol=VariableNumber)
    
    pred_mvpart <- utils::getS3method('predict', 'rpart', envir = asNamespace("mvpart"))
    
    for (i in 1:nrow(newdata)){
      if(all(is.na(newdata[i, object$Predictors]))){
        TreePrediction[i, ] <- NA
      }else{
        TreePrediction[i, ] <- pred_mvpart(object$Tree, newdata[c(i,i), ])[1, ] # There seems to be a bug that predict in "mvpart" cannot predict a single value
      }
    }
    
    colnames(TreePrediction) <- object$ResponseVariables
    
    return(TreePrediction)
    
    
  } else{
    NextMethod("predict")
  }
}


