#' Make predictions based on a multiREEMtree object on a new dataset
#'
#' Make predictions based on a multiREEMtree object on a new dataset.
#'
#' @param object A multiREEMtree object.
#' @param newdata A data frame in which to look for variables with which to predict.
#' All variables used in the fixed and random effects models, as well as the grouping
#' factors, must be included in the data frame.
#' @param level  An optional integer giving the level of grouping to be used
#' in obtaining the predictions. Defaults to be the innermost level (level=-1). If the lme.method
#' of the object is "MCMCglmm", predicting levels can only be 0 (for fixed) or 1 (for fixed + random).
#' @param ... Won't be used.
#' @return A matrix of predicted values.
#'
#' @import mvpart
#' @import MASS
#' @import expm
#' @import utils
#' @importFrom stats predict sd cov
#' @examples
#' require(multiREEMtree)
#' # simulate multivariate grouped data
#' X1 <- rnorm(100, 0, 1)
#' X2 <- rnorm(100, 0, 1)
#' X3 <- rnorm(100, 0, 1)
#' group <- rep(1:10, 10)
#' simb <- MASS::mvrnorm(10, rep(0, 2), diag(2))
#' sime <- MASS::mvrnorm(100, rep(0, 2), diag(2))
#' Y1 <- 10 * (X1>0) + 20 * (X1<0 & X2 >0) + simb[group, 1] + sime[ ,1]
#' Y2 <- 30 * (X1>0) + 50 * (X1<0 & X2 >0)+ simb[group, 2] + sime[ ,2]
#' data <- data.frame(cbind(X1, X2, X3, Y1, Y2, group))
#' f <- Y1 + Y2 ~ X1 + X2 + X3
#' r <-  ~ 1 | group
#' result <- multiREEMtree(formula=f, data=data[1:90, ], random=r)
#' pred <- predict(result, newdata=data[91:100, ], level=1)
#' @rdname predict.multiREEMtree
#' @export

predict.multiREEMtree <- function(object, newdata, level=-1, ...){
   if(class(object)=="multiREEMtree"){

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

     if(level==0) return(TreePrediction)

     if(level==-1) level = length(object$RandomEffects)

     RandomPrediction <- TreePrediction

     GV <- unlist(strsplit(object$GroupVariables, split = "\\/"))

     if(object$lme.algorithm=="MCMCglmm" & length(GV) > 1)
       stop("Predicting level of multiREEMtree using MCMCglmm can only be 0 (for fixed tree) or 1 (for fixed tree + random effect).")
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

     return(RandomPrediction)


   } else{
     NextMethod("predict")
   }
}


