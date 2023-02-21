#' Generalized Multivariate Random Effects EM Tree for longitudinal and clustered data
#'
#' \code{GmultiREEMtree} This program implements the generalized multivariate RE-EM tree method in Jing and Simonoff (2022) for 
#' non-Gaussian response variables. It is assumed that E(Y|X, b)=h(eta) and eta=f(X)+Zb, 
#' where Y is the response variable, h is a given function, f(X) is a multivariate regression tree, and b is the random effect.
#'
#' @param formula A formula of the form \code{y1 + ... + ym ~ x1 + ... + xn}, where \code{y1,...,ym} are response
#' variables and \code{x1,...,xn} are predictors.
#' @param data A data frame that contains all of the variables named in the formula.
#' @param random  A one-sided formula of the form \code{~ z1 + ... + zn | g1/.../gm} with \code{z1 + ... + zn} specifying the model
#' for the random effects and \code{g1/.../gm} the grouping structure (\code{m} may be equal to 1, in which case no / is required).
#' \code{z1,...,zn} and \code{g1,...,gm} should be contained in data. Use \code{~ 1 | g1/.../gm} to model only a random intercept for each group
#' and use \code{~ 0 + z1 + ... + zn | g1/.../gm} to remove the random intercept. 
#' @param subset  A logical vector showing the subset of data to use for fitting. Defaults to be a vector of TRUEs.
#' @param link A link function connnecting the expectation of the response variable and a linear predictor. Can either be \code{"logit"}, 
#' \code{"log"} or a self-defined list containing three elements: \code{"h"} -- a univariate function that satisfies E(Y|X, b)=h(eta); 
#' \code{"dh"} -- the first-order derivative of h; \code{"h_inv"} -- the inverse of h.
#' @param InitialEta A vector of initial estimations for eta. Defaults to be "default," which sets eta=0.75 if y=1 and eta=0.25 if y=0.
#' @param InitialRandomEffects A vector of initial values for random effects, with default value 0.
#' @param ErrorTolerance The error tolerance for the stopping criterion of the algorithm. When the
#' increase of loglikelihood for "nlme" is less than \code{ErrorTolerance},
#' the algorithm will stop. Defaults to be 0.001.
#' @param MaxIterations The largest number of iterations. Defaults to be 1000.
#' @param tree.xv Selection of tree by cross-validation: \code{"1se"} (default) - gives the simplest tree that has cross-validation error within
#' one SE of the minimum value using k-fold cross-validation, with k given by tree.xval; \code{"min"} - the tree with minimum error
#' of k-fold cross-validation, with k given by tree.xval; \code{"none"} - no cross-validation.
#' @param tree.xval Number of cross-validation folds or vector defining cross-validation groups. Defaults to be 10.
#' Could also be \code{"LOOCV"}, which will use leave-one-out cross-validation.
#' @param tree.control A list of control values used to control the \link[mvpart]{rpart}; See \link[mvpart]{rpart.control}
#' for details. Defaults to be \code{rpart.control(cp=0.01, maxsurrogate=5, usesurrogate=2)}.
#' @param autoCorrelation An optional string specifying the within-object correlation structure. Should be either NULL (default) or the name
#' of one of the available corStruct classes in \link[nlme]{corClasses} (e.g., "corAR1"). This auto-correlation structure will be applied to all of
#' the objects and response variables. 
#' @param lme.method The method used in \link[nlme]{lme}. Defaults to be \code{"REML"}. Can also be \code{"ML"}.
#' @param  lme.na.action A function that indicates what should happen when the data contain NAs. The default action (na.fail) causes lme
#' to print an error message and terminate if there are any incomplete observations. See \link[nlme]{lme}
#' for other options.
#' @param lme.control A list of control values used to control the \link[nlme]{lme} algorithm, with default opt="optim".
#' See \link[nlme]{lmeControl} for details. 
#'
#' @return A GmultiREEMtree object, which includes
#' \itemize{
#'  \item Tree - The fitted multivariate tree.
#'  \item EffectModel - The fitted linear mixed effects model.
#'  \item RandomEffects - A list of the random effects of each response variable and
#'        each random level.
#'  \item fitted - A list of the fitted values for the expectation of the responses
#'        of, with all of the fitted levels. The j-th item is a matrix equivalent to the output of using
#'        \code{fitted()} with level=j-1. (level=-1 represents the innermost level.)
#'  \item ResponseVariables - The response variables \code{y1,...,ym} defined by \code{formula}.
#'  \item RandomPredictors - The random predictors \code{z1,...,zn} defined by \code{random}.
#'  \item GroupVariables - Group Variables \code{g1/.../gm} in \code{random}.
#'  \item data - The data set.
#'  \item logLik - The final loglikelihood.
#'  \item IterationsUsed - The number of iterations used in the algorithm.
#'  \item ErrorTolerance - The parameter \code{ErrorTolerance} in the function.
#'  \item link - The parameter \code{link} in the function.
#'  \item Formula - The parameter \code{formula} in the function.
#'  \item Random - The parameter \code{random} in the function.
#'  \item Subset - The parameter \code{Subset} in the function.
#'  \item tree.control- \code{tree.control} in the function.
#'  \item xv - The parameter \code{tree.xv} in the function.
#'  \item xval - The parameter \code{tree.xval} in the function.
#'  \item lme.method - The parameter \code{lme.method} in the function.
#'  \item autoCorrelation - The fitted autoCorrelation structure.
#'  \item lme.control- The parameter \code{lme.control} in the function.
#'  \item lme.na.action - The parameter \code{lme.na.action} in the function.
#'  \item Y_original - The original values of the response variables.
#' }
#'
#' @references
#'
#' Jing, W. and Simonoff, J. S.  (2022). "A Regression Tree Method for Longitudinal and Clustered Data With Multivariate Responses." arXiv.
#'
#' Sela, R. J. and Simonoff, J. S. (2012). "RE-EM Trees: A Data Mining Approach for Longitudinal and Clustered Data." *Machine Learning*, **86**, 169-207.
#'
#' @import MASS
#' @import mvpart
#' @import nlme
#' @import Formula
#' @import expm
#' @importFrom stats as.formula logLik residuals terms na.fail cov sd
#'
#' @examples
#' require(multiREEMtree)
#' # simulate multivariate grouped data with binary response
#' X1 <- rnorm(1000, 0, 1)
#' X2 <- rnorm(1000, 0, 1)
#' X3 <- rnorm(1000, 0, 1)
#' group <- rep(1:10, 100)
#' simb <- MASS::mvrnorm(10, rep(0, 2), 0.01 * diag(2))
#' sime <- MASS::mvrnorm(1000, rep(0, 2), 0.001 * diag(2))
#' eta1 <- 2 * (X1>0) + (-2) * (X1<0 & X2 >0) + simb[group, 1] + sime[ ,1]
#' eta2 <- 2 * (X1>0) + (-2) * (X1<0 & X2 >0) + simb[group, 2] + sime[ ,2]
#' Y1 <- rep(NA, 1000)
#' Y2 <- rep(NA, 1000)
#' mu1 <- exp(eta1) / (exp(eta1) + 1)
#' mu2 <- exp(eta2) / (exp(eta2) + 1)
#' for (i in 1:1000) {
#'   Y1[i] <- rbinom(1, 1, mu1[i])
#'   Y2[i] <- rbinom(1, 1, mu2[i])
#' }
#' data <- data.frame(cbind(X1, X2, X3, Y1, Y2, group))
#' f <- Y1 + Y2 ~ X1 + X2 + X3
#' r <-  ~ 1 | group
#' result <- GmultiREEMtree(formula=f, data=data, random=r)
#' @export

GmultiREEMtree <- function(formula,
                          data,
                          random,
                          subset=NULL,
                          link="logit",
                          InitialEta="default",
                          InitialRandomEffects=0,
                          ErrorTolerance=0.001,
                          MaxIterations=1000,
                          tree.xv="1se",
                          tree.xval=10,
                          tree.control=rpart.control(cp=0.01, maxsurrogate=5, usesurrogate=2),
                          autoCorrelation=NULL,
                          lme.method="REML",
                          lme.na.action=na.fail,
                          lme.control=lmeControl(opt = "optim")
                          ){

  ## Set the initial Random Effects (default to be 0)

  if(InitialRandomEffects==0){
    InitialRandomEffects <- rep(0, dim(data)[1])
  }

  ## Subset the data if necessary

  if(identical(subset, NULL)){
    subs <- rep(TRUE, dim(data)[1])
  } else {
    if(!is.logical(subset)) stop("Subset should be a logical vector.")
    if(length(subset)!=dim(data)[1]) stop("The length of subset should be identical to the number of the observations.")
    subs <- subset
  }


  ## Parse the formula and extract the name, number and values of response variables
  if(!is.Formula(formula))
    formula <- Formula(formula)
  ResponseVariableFrame <- model.part(formula, data = data, lhs = 1)
  ResponseVariables <- colnames(ResponseVariableFrame)
  VariableNumber <- length(ResponseVariables)


  if(VariableNumber == 1) stop("For data with single response, please use REEMtree instead.")

  # Record the original response values

  OriginalResponse <- as.matrix(ResponseVariableFrame)
  
  # Define the link functions
  
  if(link=="logit") {
    h <- function(x) exp(x) / (exp(x) + 1)
    dh <- function(x) exp(x) / (exp(x) + 1) ^ 2
    h_inv <- function(x) log(x / (1-x))
    link <- list("h"=h, "dh"=dh, "h_inv"=h_inv)
  }else if(link=="log"){
    h <- function(x) exp(x) 
    dh <- function(x) exp(x) 
    h_inv <- function(x) log(x)
    link <- list("h"=h, "dh"=dh, "h_inv"=h_inv)
  }else{
    h <- link$h
    dh <- link$dh
    h_inv <- link$h_inv
  }
    
    
  # Define the initial estimator of eta
  
  if(InitialEta=="default"){
    InitialEta <- h_inv(0.25 * (OriginalResponse == 0) + 0.75 * (OriginalResponse == 1))
  }
  
  # Compute the initial pseudo responses

  PseudoResponse <- InitialEta + (1 / dh(InitialEta)) * (OriginalResponse - h(InitialEta))

  ## Extract the Predictors and the GroupVariables from the formula.

  Predictors <- setdiff(attr(terms(formula),"term.labels"), ResponseVariables)

   # The GroupVariables may be nested by /. Should split them.

  GroupVariables <- as.character(random[[2]])[[3]]
  GroupVariablesSplitted <- unlist(strsplit(GroupVariables , split = "/"))

  ## Convert the random formula into a proper form for the multivariate linear mixed effects model

    # remove all the white spaces in the formula

  RandomPredictorsFormula <- gsub(" ", "", as.character(random[[2]])[[2]], fixed = TRUE)
  RandomPredictors <- unlist(strsplit(RandomPredictorsFormula, split = "\\+"))

  if(all(RandomPredictors == "1")){

      RandomFormula <- as.formula(paste("~ 0 + trait |", GroupVariables))

      # Record the name of the Rnadom Effect terms that will be used in the output

      RandomEffectName <- ResponseVariables

  }else if("0" %in% RandomPredictors){

    # The random intercept will be canceled in this case

      RandomPredictors <- setdiff(RandomPredictors, "0")
      RandomFormula <- as.formula(paste("~ 0 + ", paste("trait:", RandomPredictors, collapse = "+"), "|", GroupVariables))
      RandomEffectName <- paste(ResponseVariables, ":", rep(RandomPredictors, each=VariableNumber), sep="")


  }else{

      RandomFormula <- as.formula(paste("~ 0 + trait + ", paste("trait:", RandomPredictors, collapse = "+"), "|", GroupVariables))
      RandomEffectName <- c(paste(ResponseVariables, ":", 1, sep=""),
                            paste(ResponseVariables, ":", rep(RandomPredictors, each=VariableNumber), sep=""))

  }


  ## Remove rows with missing response, missing group or missing all of the predictors

  missing_rows <- apply(data[ ,c(ResponseVariables, GroupVariablesSplitted)], 1, function(x){any(is.na(x))}) |
                  apply(data[,c(Predictors)], 1,  function(x){all(is.na(x))})

  if(any(missing_rows)){
    warning(paste("Row ",paste(which(missing_rows), collapse=',')," has been removed due to missing response, group or missing all of the predictors", sep = ""))
  }


  ## Update the Subset Vector including removing the missing rows

  SubsetVector <- subs & !missing_rows


  ## Set the number of CV to be the same as the length of SubsetVector If the CV method is Leave-one-out-CV


  if(tree.xval=="LOOCV") tree.xval <- sum(SubsetVector)

  tree.control$xval <- tree.xval


  ## Initial values for iteration

  ContinueCondition <- TRUE

  iterations <- 0

  AdjustedResponse <- (PseudoResponse - InitialRandomEffects)[SubsetVector, ]

  newdata <- data[SubsetVector, ]

  AdjustedRV <- paste("Adjusted.", ResponseVariables, sep="")

  oldlik <- -Inf

  trait <- as.factor(rep(1:VariableNumber, each=nrow(newdata)))

  while(ContinueCondition){
      
      value <- c(PseudoResponse)

      newdata[  ,AdjustedRV] <- AdjustedResponse

      iterations <- iterations + 1

      # Compute current tree

      tree <- mvpart::mvpart(form=formula(paste(c("data.matrix(newdata[ , AdjustedRV])", paste(Predictors, collapse = "+")),
                                          collapse = "~")),
                            data=newdata,
                            control=tree.control,
                            xv=tree.xv,
                            xval=tree.xval,
                            plot.add=F,
                            text.add=F)

      newdata[ , "nodeInd"] <- as.factor(tree$where)

      ## Estimate New Random Effects and Errors using LME

        # Get variables that identify the node for each observation

      newdataStacked <- cbind(value, newdata[rep(1:nrow(newdata), VariableNumber), ], trait)


        # Convert the auto-correlation structure

      lme.correlation <- NULL

      if(!is.null(autoCorrelation)){

        commands = paste('lme.correlation =', autoCorrelation, '(form = ~ 1 |', GroupVariables, '/trait)')

        eval(parse(text = commands))
      }

        # Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)


      if(min(tree$where)==max(tree$where)){

        # Check that the fitted tree has at least two nodes.

        lmefit <- lme(value ~ 0 + trait, data=newdataStacked, random=RandomFormula,
                      method=lme.method, control=lme.control,
                      correlation=lme.correlation, na.action=lme.na.action)

      } else {

        lmefit <- lme(value ~ 0 + trait:nodeInd, data=newdataStacked, random = RandomFormula,
                      method=lme.method, control=lme.control,
                      correlation=lme.correlation, na.action=lme.na.action)

      }


      # Get the likelihood to check on convergence
      newlik <- logLik(lmefit)
      ContinueCondition <- (newlik-oldlik > ErrorTolerance & iterations < MaxIterations)
      oldlik <- newlik

      # Update the adjusted responses
      AllEffects <- lmefit$residuals[ ,1] - lmefit$residuals[ ,dim(lmefit$residuals)[2]]
      
      
      ## Replace the predicted value on each leaf with the fixed effect coefficient in the linear mixed effects model
      adjtarg <- unique(cbind(tree$where, lmefit[["fitted"]][ ,'fixed'], trait))
      
      tree$frame$yval2[adjtarg[adjtarg[ ,3]==1, 1], ] <- matrix(adjtarg[ , 2], ncol=VariableNumber)
      
      ## Extract the fitted values. 
      
      level <- ncol(lmefit$fitted)
      FittedMatrix <- matrix(NA, ncol=VariableNumber, nrow=nrow(data), dimnames = list(1:nrow(data), ResponseVariables))
      FittedMatrix[SubsetVector, ] <- matrix(as.numeric(lmefit$fitted[ , level]),
                                              ncol=VariableNumber)

      
      eta <- FittedMatrix
      
      PseudoResponse <- eta + (1 / dh(eta)) * (OriginalResponse - h(eta))

      AdjustedResponse <- PseudoResponse - matrix(as.numeric(AllEffects), ncol=VariableNumber)

    }

  ## Replace the predicted value on each leaf with the fixed effect coefficient in the linear mixed effects model
  adjtarg <- unique(cbind(tree$where, lmefit[["fitted"]][ ,'fixed'], trait))

  tree$frame$yval2[adjtarg[adjtarg[ ,3]==1, 1], ] <- matrix(adjtarg[ , 2], ncol=VariableNumber)

  ## Extract the random effects. 
  RandomEffects <- ranef(lmefit)
  LevelNumber <- ncol(lmefit$residuals) - 1

  if(LevelNumber == 1){
      colnames(RandomEffects) <- RandomEffectName
      RandomEffects <- list(RandomEffects)
      names(RandomEffects) <- GroupVariables
  }else{
      for (l in 1:LevelNumber) {
        colnames(RandomEffects[[l]]) <- RandomEffectName
      }
  }

    ## Extract the fitted probabilities.
    
  fitted <- list()

  for(col in 1:(LevelNumber+1)){
      FittedMatrix <- matrix(NA, ncol=VariableNumber, nrow=nrow(data), dimnames = list(1:nrow(data), ResponseVariables))
      FittedMatrix[SubsetVector, ] <- matrix(as.numeric(lmefit$fitted[ , col]),
                                                ncol=VariableNumber)
      fitted_probs <- h(FittedMatrix)
      fitted <- c(fitted, list(fitted_probs))
  }
    

  names(fitted) <- colnames(lmefit$fitted)


  ## Record the autocorrelation

  lme.correlation <- lmefit$modelStruct$corStruct
    
  result <- list(Tree=tree,
                   EffectModel=lmefit,
                   RandomEffects=RandomEffects,
                   Formula=formula,
                   ResponseVariables = ResponseVariables,
                   Predictors = Predictors,
                   RandomPredictors = RandomPredictors,
                   GroupVariables = GroupVariables,
                   data=data,
                   logLik=newlik,
                   IterationsUsed=iterations,
                   Random=random,
                   Subset=SubsetVector,
                   ErrorTolerance=ErrorTolerance,
                   link=link,
                   fitted=fitted,
                   lme.method=lme.method,
                   lme.control=lme.control,
                   xv=tree.xv,
                   xval=tree.xval,
                   tree.control=tree.control,
                   autoCorrelation=lme.correlation,
                   lme.na.action=lme.na.action,
                   Y_original = OriginalResponse
                   )


  class(result) <- "GmultiREEMtree"

  return(result)
}

