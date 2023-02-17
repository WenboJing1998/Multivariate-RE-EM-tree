#' Multivariate Random Effects EM Tree for longitudinal and clustered data
#'
#' \code{multiREEMtree} This program implements the multivariate RE-EM method of Jing and Simonoff (2022),
#' which extends the RE-EM algorithm of Sela and Simonoff (2012) to a multivariate set of response variables.
#'
#' @param formula A formula of the form \code{y1 + ... + ym ~ x1 + ... + xn}, where \code{y1,...,ym} are response
#' variables and \code{x1,...,xn} are predictors.
#' @param data A data frame that contains all of the variables named in the formula.
#' @param random  A one-sided formula of the form \code{~ z1 + ... + zn | g1/.../gm} with \code{z1 + ... + zn} specifying the model
#' for the random effects and \code{g1/.../gm} the grouping structure (\code{m} may be equal to 1, in which case no / is required).
#' \code{z1,...,zn} and \code{g1,...,gm} should be contained in data. Use \code{~ 1 | g1/.../gm} to model only a random intercept for each group
#' and use \code{~ 0 + z1 + ... + zn | g1/.../gm} to remove the random intercept. If \code{lme.algorithm="MCMCglmm"}, only \code{m}=1 is allowed.
#' @param subset  A logical vector showing the subset of data to use for fitting. Defaults to be a vector of TRUEs.
#' @param stand.y An option to standardize the response variables, should be one of "marginal" (standardize each response variable by subtracting
#' its mean and dividing by its sample standard deviation), "covariance"(standardize the set of response variables using the vector of sample
#' means and sample covariance matrix, i.e., convert y to y * (cov(y)^(-1/2)), or "none" (no standardization). Defaults to be "marginal".
#' @param InitialRandomEffects A vector of initial values for random effects, with default value 0.
#' @param ErrorTolerance The error tolerance for the stopping criterion of the algorithm. When the
#' increase of loglikelihood for "nlme" (or the decrease of DIC for "MCMCglmm") is less than \code{ErrorTolerance},
#' the algorithm will stop. Defaults to be 0.001.
#' @param MaxIterations The largest number of iterations. Defaults to be 1000.
#' @param tree.xv Selection of tree by cross-validation: \code{"1se"} (default) - gives the simplest tree that has cross-validation error within
#' one SE of the minimum value using k-fold cross-validation, with k given by tree.xval; \code{"min"} - the tree with minimum error
#' of k-fold cross-validation, with k given by tree.xval; \code{"none"} - no cross-validation.
#' @param tree.xval Number of cross-validation folds or vector defining cross-validation groups. Defaults to be 10.
#' Could also be \code{"LOOCV"}, which will use leave-one-out cross-validation.
#' @param tree.control A list of control values used to control the \link[mvpart]{rpart}; See \link[mvpart]{rpart.control}
#' for details. Defaults to be \code{rpart.control(cp=0.01, maxsurrogate=5, usesurrogate=2)}.
#' @param lme.algorithm The algorithm used to fit linear mixed effect model. Either \code{"nlme"} or \code{"MCMCglmm"}. Defaults to be "nlme".
#' @param autoCorrelation An optional string specifying the within-object correlation structure. Should be either NULL (default) or the name
#' of one of the available corStruct classes in \link[nlme]{corClasses} (e.g., "corAR1"). This auto-correlation structure will be applied to all of
#' the objects and response variables. Only used when \code{lme.algorithm="nlme"}.
#' @param lme.method The method used in \link[nlme]{lme}. Defaults to be \code{"REML"}. Can also be \code{"ML"}. Only used when \code{lme.algorithm="nlme"}.
#' @param  lme.na.action A function that indicates what should happen when the data contain NAs. The default action (na.fail) causes lme
#' to print an error message and terminate if there are any incomplete observations. Only used when \code{lme.algorithm="nlme"}. See \link[nlme]{lme}
#' for other options.
#' @param lme.control A list of control values used to control the \link[nlme]{lme} algorithm, with default opt="optim".
#' See \link[nlme]{lmeControl} for details. Only used when \code{lme.algorithm="nlme"}.
#'
#' @return A multiREEMtree object, which includes.
#' \itemize{
#'  \item Tree - The fitted multivariate tree.
#'  \item EffectModel - The fitted linear mixed effects model.
#'  \item RandomEffects - A list of the random effects of each response variable and
#'        each random level.
#'  \item fitted - A list of the fitted values of all of the levels.
#'        The j-th item is a matrix equivalent to the output of using
#'        \code{fitted()} with level=j-1. (level=-1 represents the innermost level.)
#'  \item residuals - A list of the residual values of all of the levels. The
#'        j-th item is a matrix equivalent to the output of using
#'        \code{residual()} with level=j-1.
#'  \item ResponseVariables - The response variables \code{y1,...,ym} defined by \code{formula}.
#'  \item RandomPredictors - The random predictors \code{z1,...,zn} defined by \code{random}.
#'  \item GroupVariables - Group Variables \code{g1/.../gm} in \code{random}.
#'  \item data - The data set.
#'  \item logLik - The final loglikelihood if \code{lme.algorithm="nlme"}.
#'  \item DIC - The final value of Deviance Information Criterion if \code{lme.algorithm="MCMCglmm"}.
#'  \item IterationsUsed - The number of iterations used in the algorithm.
#'  \item ErrorTolerance - The parameter \code{ErrorTolerance} in the function.
#'  \item Formula - The parameter \code{formula} in the function.
#'  \item Random - The parameter \code{random} in the function.
#'  \item Subset - The parameter \code{Subset} in the function.
#'  \item tree.control- \code{tree.control} in the function.
#'  \item xv - The parameter \code{tree.xv} in the function.
#'  \item xval - The parameter \code{tree.xval} in the function.
#'  \item lme.algorithm - The parameter \code{lme.algorithm} in the function.
#'  \item lme.method - The parameter \code{lme.method} in the function.
#'  \item autoCorrelation - The fitted autoCorrelation structure.
#'  \item lme.control- The parameter \code{lme.control} in the function.
#'  \item lme.na.action - The parameter \code{lme.na.action} in the function.
#'  \item Y_original - The original values of the response variables.
#'  \item stand.y - The parameter \code{stand.y} in the function.
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
#' @import MCMCglmm
#' @import expm
#' @importFrom stats as.formula logLik residuals terms na.fail cov sd
#'
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
#' result <- multiREEMtree(formula=f, data=data, random=r)
#' @export

multiREEMtree <- function(formula,
                          data,
                          random,
                          subset=NULL,
                          stand.y = "marginal",
                          InitialRandomEffects=0,
                          ErrorTolerance=0.001,
                          MaxIterations=1000,
                          tree.xv="1se",
                          tree.xval=10,
                          tree.control=rpart.control(cp=0.01, maxsurrogate=5, usesurrogate=2),
                          lme.algorithm="nlme",
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

  ## Record the original response values

  OriginalResponse <- as.matrix(ResponseVariableFrame)



  ## Standardize the responses

  if(! stand.y %in% c("covariance", "marginal","none")) stop("stand.y should be covariance, marginal or none.")

  if(stand.y == "covariance"){

    cov.y <- cov(ResponseVariableFrame)
    StandardizedResponse <- OriginalResponse %*% sqrtm(solve(cov.y))

    # A function to convert y from the standardized version to the original version

    standard_to_original <- function(y) {
      y %*% sqrtm(cov.y)
    }

    # A function to convert the random effect from the standardized version to the original version

    standard_to_original_RE <- function(RE) {
      as.matrix(RE) %*% sqrtm(cov.y)
    }

  }else if(stand.y == "marginal"){

    y_mean <- apply(OriginalResponse, 2, mean)

    y_std <- apply(OriginalResponse, 2, sd)

    StandardizedResponse <- t((t(OriginalResponse) - y_mean) / y_std)

    standard_to_original <- function(y) {
      t(t(y) * y_std + y_mean)
    }

    standard_to_original_RE <- function(RE) {
      t(t(as.matrix(RE) * y_std))
    }

  }else{

    StandardizedResponse <- OriginalResponse

    standard_to_original <- function(y) y

    standard_to_original_RE <- function(RE) RE


  }

  ## Replace the Response Variables with the standardized version

  data[ , ResponseVariables] <- StandardizedResponse

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

  AdjustedResponse <- (StandardizedResponse - InitialRandomEffects)[SubsetVector, ]

  newdata <- data[SubsetVector, ]

  AdjustedRV <- paste("Adjusted.", ResponseVariables, sep="")

  if(lme.algorithm=="nlme"){

    oldlik <- -Inf

    # Make a new data frame to include all the new variables

    value <- c(as.matrix(newdata[ ,ResponseVariables]))
    trait <- as.factor(rep(1:VariableNumber, each=nrow(newdata)))


    while(ContinueCondition){


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

      AdjustedResponse <- newdata[ ,ResponseVariables] - matrix(as.numeric(AllEffects), ncol=VariableNumber)

    }

    ## Replace the predicted value on each leaf with the fixed effect coefficient in the linear mixed effects model
    adjtarg <- unique(cbind(tree$where, lmefit[["fitted"]][ ,'fixed'], trait))

    tree$frame$yval2[adjtarg[adjtarg[ ,3]==1, 1], ] <- standard_to_original(matrix(adjtarg[ , 2], ncol=VariableNumber))

    ## Extract the random effects. Use standard_to_original_RE to recover the original values.
    RandomEffects <- ranef(lmefit)
    LevelNumber <- ncol(lmefit$residuals) - 1

    if(LevelNumber == 1){
      RandomEffects <- standard_to_original_RE(RandomEffects)
      colnames(RandomEffects) <- RandomEffectName
      RandomEffects <- list(RandomEffects)
      names(RandomEffects) <- GroupVariables
    }else{
      for (l in 1:LevelNumber) {
        RandomEffects[[l]] <- standard_to_original_RE(RandomEffects[[l]])
        colnames(RandomEffects[[l]]) <- RandomEffectName
      }
    }

    ## Extract the fitted values. Use standard_to_original to recover the original values.

    fitted <- list()

    for(col in 1:(LevelNumber+1)){
      FittedMatrix <- matrix(NA, ncol=VariableNumber, nrow=nrow(data), dimnames = list(1:nrow(data), ResponseVariables))
      FittedMatrix[SubsetVector, ] <- standard_to_original(matrix(as.numeric(lmefit$fitted[ , col]),
                                                ncol=VariableNumber))
      fitted <- c(fitted, list(FittedMatrix))

    }

    names(fitted) <- colnames(lmefit$fitted)


    ## Extract the residuals.

    residuals <- list()

    for(col in 1:(LevelNumber+1)){
      ResidualMatrix <- OriginalResponse -  fitted[[col]]
      residuals <- c(residuals, list(ResidualMatrix))

    }

    names(residuals) <- colnames(lmefit$residuals)

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
                   lme.algorithm="nlme",
                   fitted=fitted,
                   residuals=residuals,
                   lme.method=lme.method,
                   lme.control=lme.control,
                   xv=tree.xv,
                   xval=tree.xval,
                   tree.control=tree.control,
                   autoCorrelation=lme.correlation,
                   lme.na.action=lme.na.action,
                   Y_original = OriginalResponse,
                   stand.y = stand.y
                   )


  }else if(lme.algorithm=="MCMCglmm"){

    if(length(GroupVariablesSplitted) > 1)

      stop("Please set lme.algorithm to be 'nlme' for multi-level randomness.")

    # Modify formulae for MCMC

    FixedFormulaMCMC <- as.formula(paste("cbind(", paste(ResponseVariables, collapse=","), ")~0 + trait:nodeInd"))

    if(all(RandomPredictors == "1")){

      RandomFormulaMCMC <- as.formula(paste("~us(trait):",GroupVariablesSplitted, sep=""))

    }else{

      RandomFormulaMCMC <- as.formula(paste(c("~us(trait):",
                                              paste("us(trait:", RandomPredictors, "):", sep="")),
                                            GroupVariablesSplitted,
                                            sep="",
                                            collapse="+"))
    }


    oldDIC <- Inf

    while(ContinueCondition){

      # Current values of variables
      newdata[ ,AdjustedRV] <- AdjustedResponse
      iterations <- iterations + 1

      # Compute current tree

      tree <- mvpart::mvpart(form=formula(paste(c("data.matrix(newdata[ ,AdjustedRV])", paste(Predictors, collapse = "+")),
                                                collapse = "~")),
                             data=newdata,
                             control=tree.control,
                             xv=tree.xv,
                             xval=tree.xval,
                             plot.add=F,
                             text.add=F)


      newdata[ , "nodeInd"] <- as.factor(tree$where)


      MCMCglmmfit <- MCMCglmm(fixed = FixedFormulaMCMC,
                              random = RandomFormulaMCMC,
                              rcov = ~ idh(trait):units,
                              pr = TRUE,
                              family = rep("gaussian", VariableNumber),
                              data = newdata,
                              verbose=FALSE)


      # Use DIC instead of loglik in MCMC models
      newDIC <- MCMCglmmfit$DIC
      ContinueCondition <- (oldDIC-newDIC > ErrorTolerance & iterations < MaxIterations)
      oldDIC <- newDIC

      # Compute the random effects of each sample fitted by MCMCglmm

      FittedFixedMCMC <- predict(MCMCglmmfit, newdata, marginal=MCMCglmmfit$Random$formula)
      FittedAllMCMC <- predict(MCMCglmmfit, newdata, marginal=NULL)


      # Update the adjusted response variables
      AdjustedResponse <- newdata[ , ResponseVariables] - (FittedAllMCMC - FittedFixedMCMC)


    }



    # Replace the predicted value on each leaf with the fixed effect coefficient in the linear mixed effects model
    adjtarg <- unique(cbind(tree$where, matrix(FittedFixedMCMC, ncol=2, byrow=F)))

    tree$frame$yval2[adjtarg[ , 1], ] <- standard_to_original(adjtarg[ , -1])


    EffectsMCMC <- as.matrix(apply(MCMCglmmfit[["Sol"]], 2, mean))

    NodesNumber <- length(unique(newdata[ ,"nodeInd"]))
    RandomEffectsMCMC <- matrix(EffectsMCMC[-(1 : (VariableNumber * NodesNumber))],
                                 byrow = F,
                                 ncol = length(RandomEffectName))
    colnames(RandomEffectsMCMC) <- RandomEffectName

    # Change the rownames of RandomEffectsMCMC into the groupnames
    GroupNumber <- length(unique(newdata[ ,GroupVariables]))
    RownamesTemp <- as.matrix(rownames(EffectsMCMC)[-(1 : (VariableNumber * NodesNumber))])
    rownames(RandomEffectsMCMC) <- apply(as.matrix(RownamesTemp[c(1:GroupNumber), ]), 1,
                                         function(x) strsplit(x, paste(".",GroupVariablesSplitted,".", sep=""))[[1]][2])


    ## Extract the Random Effects, fitted values and residuals

    RandomEffectsMCMC <- standard_to_original_RE(RandomEffectsMCMC)

    RandomEffectsMCMC <- list(RandomEffectsMCMC)
    names(RandomEffectsMCMC) <- GroupVariablesSplitted



    FittedMatrixFixedMCMC <- matrix(NA, nrow=nrow(data), ncol=VariableNumber, dimnames = list(1:nrow(data), ResponseVariables))

    FittedMatrixFixedMCMC[SubsetVector, ] <- standard_to_original(matrix(FittedFixedMCMC, ncol=VariableNumber))

    FittedMatrixRandMCMC <- matrix(NA, nrow=nrow(data), ncol=VariableNumber, dimnames = list(1:nrow(data), ResponseVariables))

    FittedMatrixRandMCMC[SubsetVector, ] <- standard_to_original((matrix(FittedAllMCMC, ncol=VariableNumber)))

    fitted <- list(fixed=FittedMatrixFixedMCMC, rand=FittedMatrixRandMCMC)

    names(fitted)[2] <- GroupVariablesSplitted

    residuals <- list()

    for(col in 1:2){
      ResidualMatrix <- OriginalResponse -  fitted[[col]]
      residuals <- c(residuals, list(ResidualMatrix))

    }


    names(residuals) <- names(fitted)


    result <- list(Tree=tree,
                   EffectModel=MCMCglmmfit,
                   RandomEffects=RandomEffectsMCMC,
                   Formula=formula,
                   Random=random,
                   ResponseVariables = ResponseVariables,
                   Predictors = Predictors,
                   RandomPredictors = RandomPredictors,
                   GroupVariables = GroupVariables,
                   data=data,
                   DIC=newDIC,
                   IterationsUsed=iterations,
                   Subset=SubsetVector,
                   ErrorTolerance=ErrorTolerance,
                   fitted = fitted,
                   residuals=residuals,
                   lme.algorithm="MCMCglmm",
                   tree.control=tree.control,
                   xv=tree.xv,
                   xval=tree.xval,
                   Y_original = OriginalResponse,
                   stand.y = stand.y)

  }else{
    stop("lme.algorithm should be either nlme or MCMCglmm.")
  }


  class(result) <- "multiREEMtree"

  return(result)
}

