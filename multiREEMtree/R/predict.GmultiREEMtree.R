#' Make predictions based on a GmultiREEMtree object on a new dataset
#'
#' Make predictions based on a GmultiREEMtree object on a new dataset.
#'
#' @param object A multiREEMtree object.
#' @param newdata A data frame in which to look for variables with which to predict.
#' All variables used in the fixed and random effects models, as well as the grouping
#' factors, must be included in the data frame.
#' @param level  An optional integer giving the level of grouping to be used
#' in obtaining the predictions. Defaults to be the innermost level (level=-1). If the lme.method
#' of the object is "MCMCglmm", predicting levels can only be 0 (for fixed) or 1 (for fixed + random).
#' @param ... Won't be used.
#' @return A matrix of predicted values for the expectation of the responses.
#'
#' @import mvpart
#' @import MASS
#' @import expm
#' @import utils
#' @importFrom stats predict sd cov
#' @rdname predict.GmultiREEMtree
#' @export
predict.GmultiREEMtree <- function(object, newdata, level=-1, ...){
   if(inherits(object, "GmultiREEMtree")){
     

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
     
     h <- object$link$h

     if(level==0) return(h(TreePrediction))

     if(level==-1) level = length(object$RandomEffects)

     RandomPrediction <- TreePrediction

     GV <- unlist(strsplit(object$GroupVariables, split = "\\/"))

     if(level > length(object$RandomEffects))
       stop(paste("The max level for this object is", length(object$RandomEffects)))

     for (i in 1:nrow(newdata)) {


       if(all(object$RandomPredictors == "1")){

         RandomPredictorsValue <- rep(1, VariableNumber)

       }else if(colnames(object$RandomEffects[[level]])[1] == paste(object$ResponseVariables[1], ":", "1", sep="")){

         # An intercept is included
         RandomPredictorsValue <- rep(as.numeric(c(1, newdata[i, object$RandomPredictors])), each = VariableNumber)

       }else{

         #No intercept
         RandomPredictorsValue <- rep(as.numeric(newdata[i, object$RandomPredictors]), each = VariableNumber)

       }

         group <- paste(newdata[i, GV][1:level], collapse="/")

         if(group  %in% rownames(object$RandomEffects[[level]])){
           rpred <- as.numeric(object$RandomEffects[[level]][group,  ] *  RandomPredictorsValue)
           RandomPrediction[i,  ] <- as.numeric(TreePrediction[i,  ] +
                                                  apply(matrix(rpred, ncol = VariableNumber, byrow=T), 2, sum))
         }

     }

     colnames(RandomPrediction) <- object$ResponseVariables
     rownames(RandomPrediction) <- rownames(newdata)

     return(h(RandomPrediction))


   } else{
     NextMethod("predict")
   }
}


