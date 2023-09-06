  # R code for the simulation studies of multivariate RE-EM tree with no random effect. 
  # The code is run on NYU High Performance Computing (NYU HPC).

  library(MASS, warn.conflicts = F, quietly = T)
  library(mvpart, warn.conflicts = F, quietly = T)
  library(MCMCglmm, warn.conflicts = F, quietly = T)
  library(nlme, warn.conflicts = F, quietly = T)
  library(foreach, warn.conflicts = F, quietly = T)
  library(doParallel, warn.conflicts = F, quietly = T)
  library(expm, warn.conflicts = F, quietly = T)
  library(Formula, warn.conflicts = F, quietly = T)
  library(doRNG, warn.conflicts = F, quietly = T)

  REEMtree <- function(formula, data, random, subset=NULL, initialRandomEffects=rep(0, TotalObs),
                       ErrorTolerance=0.001, MaxIterations=1000, verbose=FALSE, tree.control=rpart.control(),
                       cv=TRUE, cpmin = 0.001, no.SE = 1,
                       lme.control=lmeControl(opt = "optim"), method="REML", correlation=NULL){
     TotalObs <- dim(data)[1]
     
     originaldata <- data
     
     # Subset the data if necessary
     if(identical(subset, NULL)){
        subs <- rep(TRUE, dim(data)[1])
     } else {
        subs <- subset
     }
     
     # Parse formula
     Predictors <- paste(attr(terms(formula),"term.labels"),collapse="+")
     TargetName <- formula[[2]]
     # Remove the name of the data frame if necessary
     if(length(TargetName)>1)
        TargetName <-TargetName[3]
     if(verbose) print(paste("Target variable: ", TargetName))
     
     Target <- data[,toString(TargetName)]
     
     
     # Condition that indicates the loop has not converged or
     # run out of iterations
     ContinueCondition <- TRUE
     
     iterations <- 0
     
     # Get initial values
     AdjustedTarget <- Target - initialRandomEffects
     oldlik <- -Inf
     
     # Make a new data frame to include all the new variables
     newdata <- data
     newdata[, "SubsetVector"] <- subs
     
     while(ContinueCondition){
        
        # Current values of variables
        newdata[,"AdjustedTarget"] <- AdjustedTarget
        iterations <- iterations+1
        
        # Compute current tree
        if (cv) {
           
           tree1 <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                        collapse = "~")), data = newdata, subset = subs,
                          method = "anova", control = rpart.control(cp=cpmin))
           if (nrow(tree1$cptable)==1){
              tree <- tree1}
           else {
              cventry <- which.min(tree1$cptable[, "xerror"])
              if (no.SE == 0){
                 cpcv <- tree1$cptable[cventry, "CP"]
                 tree <- prune(tree1, cp=cpcv)}
              else {
                 xerrorcv <- tree1$cptable[cventry, "xerror"]
                 sexerrorcv <- xerrorcv + tree1$cptable[cventry, "xstd"] * no.SE
                 cpcvse <- tree1$cptable[which.max(tree1$cptable[, "xerror"] <= sexerrorcv), "CP"]
                 tree <- prune(tree1, cp=cpcvse)}
           }
        }
        else {
           tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
                                       collapse = "~")), data = newdata, subset = subs,
                         method = "anova", control = rpart.control())
           
           # tree <- rpart(formula(paste(c("AdjustedTarget", Predictors),
           #      collapse = "~")), data = newdata, subset = subs,
           #      method = "anova", control = tree.control)
        }
        if(verbose) print(tree)
        
        ## Estimate New Random Effects and Errors using LME
        # Get variables that identify the node for each observation
        newdata[ ,"nodeInd"] <- 0
        newdata[subs,"nodeInd"] <- tree$where
        # Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)
        # Check that the fitted tree has at least two nodes.
        if(min(tree$where)==max(tree$where)){
           lmefit <- lme(formula(paste(c(toString(TargetName),1), collapse="~")), data=newdata, random=random,
                         subset=SubsetVector, method=method, control=lme.control, correlation=correlation)
        } else {
           lmefit <- lme(formula(paste(c(toString(TargetName),"as.factor(nodeInd)"), collapse="~")), data=newdata, random=random,
                         subset=SubsetVector, method=method, control=lme.control, correlation=correlation)
        }
        
        # population prediction for each leaf
        adjtarg <- unique(cbind(tree$where, predict(lmefit, level=0)))
        tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]
        
        if(verbose){
           print(lmefit)
           print(paste("Estimated Error Variance = ", lmefit$sigma))
           print("Estimated Random Effects Variance = ")
           print(as.matrix(lmefit$modelStruct$reStruct[[1]])*lmefit$sigma^2)
        }
        
        # Get the likelihood to check on convergence
        newlik <- logLik(lmefit)
        if(verbose) print(paste("Log likelihood: ", newlik))
        
        ContinueCondition <- (newlik-oldlik>ErrorTolerance & iterations < MaxIterations)
        oldlik <- newlik
        
        # Extract random effects to make the new adjusted target
        AllEffects <- lmefit$residuals[,1]-lmefit$residuals[,dim(lmefit$residuals)[2]]
        AdjustedTarget[subs] <- Target[subs] - AllEffects
     }
     
     residuals <- rep(NA, length=length(Target))
     residuals[subs] <- Target[subs]-predict(lmefit)
     attr(residuals, "label") <- NULL
     
     
     adjtarg <- unique(cbind(tree$where, predict(lmefit, level=0)))
     tree$frame[adjtarg[,1],]$yval <- adjtarg[,2]
     
     
     
     result <- list(Tree=tree, EffectModel=lmefit, RandomEffects=ranef(lmefit),
                    BetweenMatrix=as.matrix(lmefit$modelStruct$reStruct[[1]])*lmefit$sigma^2,
                    ErrorVariance=lmefit$sigma^2, data=data, logLik=newlik,
                    IterationsUsed=iterations, Formula=formula, Random=random, Subset=subs,
                    ErrorTolerance=ErrorTolerance, correlation=correlation,
                    residuals=residuals, method=method, cv=cv, lme.control=lme.control, tree.control=tree.control)
     class(result) <- "REEMtree"
     
     return(result)
  }
  
  
  multitree <- function (formula, data, subset = NULL, stand.y = "marginal", 
                         tree.xv = "1se", tree.xval = 10, tree.control = rpart.control(cp = 0.01, 
                                                                                       maxsurrogate = 5, usesurrogate = 2)) 
    #
    # Wrapper function for calling of mvpart function
    #
  {
    if (identical(subset, NULL)) {
      subs <- rep(TRUE, dim(data)[1])
    }
    else {
      if (!is.logical(subset)) 
        stop("Subset should be a logical vector.")
      if (length(subset) != dim(data)[1]) 
        stop("The length of subset should be identical to the number of the observations.")
      subs <- subset
    }
    if (!is.Formula(formula)) 
      formula <- Formula(formula)
    ResponseVariableFrame <- model.part(formula, data = data, 
                                        lhs = 1)
    ResponseVariables <- colnames(ResponseVariableFrame)
    VariableNumber <- length(ResponseVariables)
    if (VariableNumber == 1) 
      stop("For data with single response, please use CART instead.")
    OriginalResponse <- as.matrix(ResponseVariableFrame)
    if (!stand.y %in% c("covariance", "marginal", 
                        "none")) 
      stop("stand.y should be covariance, marginal or none.")
    if (stand.y == "covariance") {
      cov.y <- cov(ResponseVariableFrame)
      StandardizedResponse <- OriginalResponse %*% sqrtm(solve(cov.y))
      standard_to_original <- function(y) {
        y %*% sqrtm(cov.y)
      }
      standard_to_original_RE <- function(RE) {
        as.matrix(RE) %*% sqrtm(cov.y)
      }
    }
    else if (stand.y == "marginal") {
      y_mean <- apply(OriginalResponse, 2, mean)
      y_std <- apply(OriginalResponse, 2, sd)
      StandardizedResponse <- t((t(OriginalResponse) - y_mean)/y_std)
      standard_to_original <- function(y) {
        t(t(y) * y_std + y_mean)
      }
      standard_to_original_RE <- function(RE) {
        t(t(as.matrix(RE) * y_std))
      }
    }
    else {
      StandardizedResponse <- OriginalResponse
      standard_to_original <- function(y) y
      standard_to_original_RE <- function(RE) RE
    }
    data[, ResponseVariables] <- StandardizedResponse
    Predictors <- setdiff(attr(terms(formula), "term.labels"), 
                          ResponseVariables)
    missing_rows <- apply(data[, ResponseVariables], 
                          1, function(x) {
                            any(is.na(x))
                          }) | apply(data[, c(Predictors)], 1, function(x) {
                            all(is.na(x))
                          })
    if (any(missing_rows)) {
      warning(paste("Row ", paste(which(missing_rows), 
                                  collapse = ","), " has been removed due to missing response or missing all of the predictors", 
                    sep = ""))
    }
    SubsetVector <- subs & !missing_rows
    if (tree.xval == "LOOCV") 
      tree.xval <- sum(SubsetVector)
    tree.control$xval <- tree.xval
    
    
    tree <- mvpart::mvpart(form = formula(paste(c("data.matrix(data[ , ResponseVariables])", 
                                                  paste(Predictors, collapse = "+")), collapse = "~")), 
                           data = data, control = tree.control, xv = tree.xv, 
                           xval = tree.xval, plot.add = F, text.add = F)
    
    
    tree$frame$yval2[tree$frame$var=="<leaf>", ] <- standard_to_original(tree$frame$yval2[tree$frame$var=="<leaf>", ])
    
    result <- list(Tree = tree, Formula = formula, ResponseVariables = ResponseVariables, 
                   Predictors = Predictors, data = data, Subset = SubsetVector, 
                   xv = tree.xv, xval = tree.xval, tree.control = tree.control, Y_original = OriginalResponse, 
                   stand.y = stand.y)
    class(result) <- "multitree"
    return(result)
  }
  
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
                            correlation=NULL,
                            lme.method="REML",
                            lme.na.action=na.fail,
                            lme.control=lmeControl(opt = "optim")
  ){
    
    if(! stand.y %in% c("covariance", "marginal","none")) stop("stand.y should be covariance, marginal or none")
    
    
    if(InitialRandomEffects==0){
      InitialRandomEffects <- rep(0, dim(data)[1])
    }
    
    # Subset the data if necessary
    
    if(identical(subset, NULL)){
      subs <- rep(TRUE, dim(data)[1])
    } else {
      if(!is.logical(subset)) stop("Subset should be a logical vector.")
      if(length(subset)!=dim(data)[1]) stop("The length of subset should be identical to the number of the observations.")
      subs <- subset
    }
    
    
    # Parse formula
    if(!is.Formula(formula))
      formula <- Formula(formula)
    ResponseVariableFrame <- model.part(formula, data = data, lhs = 1)
    ResponseVariables <- colnames(ResponseVariableFrame)
    VariableNumber <- length(ResponseVariables)
    
    #PredictorsFrame <- model.part(formula, data = data, rhs = 1)
    #Predictors <- colnames(PredictorsFrame)
    
    if(VariableNumber == 1) stop("For data with single response, please use REEMtree instead.")
    
    # Standardize the responses
    
    OriginalResponse <- as.matrix(ResponseVariableFrame)
    
    
    if(stand.y == "covariance"){
      
      cov.y <- cov(ResponseVariableFrame)
      StandardizedResponse <- OriginalResponse %*% sqrtm(solve(cov.y))
      
      # Function to convert y from the standardized version to the original version
      
      standard_to_original <- function(y) {
        y %*% sqrtm(cov.y)
      }
      
      # Function to convert the random effect from the standardized version to the original version
      
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
    
    
    data[ , ResponseVariables] <- StandardizedResponse
    
    # Read the Predictors and the GroupVariables from the formula
    
    Predictors <- setdiff(attr(terms(formula),"term.labels"), ResponseVariables)
    
    GroupVariables <- as.character(random[[2]])[[3]]
    GroupVariablesSplitted <- unlist(strsplit(GroupVariables , split = "/"))
    
    # Convert the random formula into a proper form for the multivariate linear mixed effects model
    RandomPredictorsFormula <- gsub(" ", "", as.character(random[[2]])[[2]], fixed = TRUE) # remove all the white spaces in the formula
    RandomPredictors <- unlist(strsplit(RandomPredictorsFormula, split = "\\+"))
    
    if(all(RandomPredictors == "1")){
      RandomFormula <- as.formula(paste("~ 0 + trait |", GroupVariables))
      RandomEffectName <- ResponseVariables # Record the name that will be used in the output
    }else if("0" %in% RandomPredictors){
      RandomPredictors <- setdiff(RandomPredictors, "0")
      RandomFormula <- as.formula(paste("~ 0 + ", paste("trait:", RandomPredictors, collapse = "+"), "|", GroupVariables))
      RandomEffectName <- paste(ResponseVariables, ":", rep(RandomPredictors, each=VariableNumber), sep="")
    }else{
      RandomFormula <- as.formula(paste("~ 0 + trait + ", paste("trait:", RandomPredictors, collapse = "+"), "|", GroupVariables))
      RandomEffectName <- c(paste(ResponseVariables, ":", 1, sep=""),
                            paste(ResponseVariables, ":", rep(RandomPredictors, each=VariableNumber), sep=""))
    }
    
    
    # Remove rows with missing response, missing group or missing all the predictors
    
    missing_rows <- apply(data[ ,c(ResponseVariables, GroupVariablesSplitted)], 1, function(x){any(is.na(x))}) |
      apply(data[,c(Predictors)], 1,  function(x){all(is.na(x))})
    if(any(missing_rows)){
      # data <- data[!missing_rows, ]
      warning(paste("Row ",paste(which(missing_rows), collapse=',')," has been removed due to missing response, group or missing all the predictors", sep = ""))
    }
    
    
    SubsetVector <- subs & !missing_rows
    
    
    if(tree.xval=="LOOCV") tree.xval <- sum(SubsetVector)
    
    tree.control$xval <- tree.xval
    
    
    ContinueCondition <- TRUE
    
    iterations <- 0
    
    AdjustedResponse <- (StandardizedResponse - InitialRandomEffects)[SubsetVector, ]
    
    newdata <- data[SubsetVector, ]
    
    #newdata[ , "SubsetVector"] <- SubsetVector
    
    AdjustedRV <- paste("Adjusted.", ResponseVariables, sep="")
    
    if(lme.algorithm=="nlme"){
      # Get initial values
      
      oldlik <- -Inf
      
      # Make a new data frame to include all the new variables
      
      value <- c(as.matrix(newdata[ ,ResponseVariables]))
      trait <- as.factor(rep(1:VariableNumber, each=nrow(newdata)))
      
      
      while(ContinueCondition){
        
        #newdata[  ,AdjustedRV] <- OriginalResponse
        
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
        
        
        # Fit linear model with nodes as predictors (we use the original target so likelihoods are comparable)
        
        
        
        # Check that the fitted tree has at least two nodes.
        if(min(tree$where)==max(tree$where)){
          lmefit <- lme(value ~ 0 + trait, data=newdataStacked, random=RandomFormula,
                        method=lme.method, control=lme.control,
                        correlation=correlation, na.action=lme.na.action)
        } else {
          lmefit <- lme(value ~ 0 + trait:nodeInd, data=newdataStacked, random = RandomFormula,
                        method=lme.method, control=lme.control,
                        correlation=correlation, na.action=lme.na.action)
        }
        

        
        
        # Get the likelihood to check on convergence
        newlik <- logLik(lmefit)
        ContinueCondition <- (newlik-oldlik > ErrorTolerance & iterations < MaxIterations)
        oldlik <- newlik
        
        # Update the adjusted responses
        AllEffects <- lmefit$residuals[ ,1] - lmefit$residuals[ ,dim(lmefit$residuals)[2]]
        
        AdjustedResponse <- newdata[ ,ResponseVariables] - matrix(as.numeric(AllEffects), ncol=VariableNumber)
        
      }
      
      # Replace the predicted value on each leaf with the fixed effect coefficient in the linear mixed effects model
      adjtarg <- unique(cbind(tree$where, lmefit[["fitted"]][ ,'fixed'], trait))
      
      tree$frame$yval2[adjtarg[adjtarg[ ,3]==1, 1], ] <- standard_to_original(matrix(adjtarg[ , 2], ncol=VariableNumber))
      
      # Extract the random effects
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
      
      fitted <- list()
      
      for(col in 1:(LevelNumber+1)){
        FittedMatrix <- matrix(NA, ncol=VariableNumber, nrow=nrow(data), dimnames = list(1:nrow(data), ResponseVariables))
        FittedMatrix[SubsetVector, ] <- standard_to_original(matrix(as.numeric(lmefit$fitted[ , col]),
                                                                    ncol=VariableNumber))
        fitted <- c(fitted, list(FittedMatrix))
        
      }
      
      names(fitted) <- colnames(lmefit$fitted)
      
      residuals <- list()
      
      for(col in 1:(LevelNumber+1)){
        ResidualMatrix <- OriginalResponse -  fitted[[col]]
        residuals <- c(residuals, list(ResidualMatrix))
        
      }
      
      names(residuals) <- colnames(lmefit$residuals)
      
      
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
                     correlation=correlation,
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
      
      
      # Get initial values
      
      oldDIC <- Inf
      
      ContinueCondition <- TRUE
      
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
                     correlation=correlation,
                     lme.na.action=lme.na.action,
                     Y_original = OriginalResponse,
                     stand.y = stand.y)
      
    }else{
      stop("lme.algorithm should be either nlme or MCMCglmm.")
    }
    
    
    class(result) <- "multiREEMtree"
    
    return(result)
  }
  
  predict.multitree <- function(object, newdata, level=-1, ...){
    if(class(object)=="multitree"){
      
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
  
  # Function: examine if two rpart objects have the same tree structure
  is.sametree.rpart <- function(object1, object2){
    if(nrow(object1$frame) != nrow(object2$frame)) return(FALSE)
    return(
      identical(as.integer(row.names(object1$frame)), as.integer(row.names(object2$frame))) &
        identical(as.vector(object1$frame$var), as.vector(object2$frame$var)))
  }

 cl <- makeCluster(6)
 registerDoParallel(cl)
 
 # No random effect!
 I_seq <- c(50, 100, 200, 400, 800) # Number of objects
 Ti_seq <- c(5, 10, 25, 50)  # Number of observations of each objects
 sigmae_seq <- c(0.5, 1, 1.5, 2) # Sequence of the variance of the noise
 structure_seq <- c(1, 2, 3) # Sequence of tree structure: 1 for bivariate simple tree, 2 for bivariate complex tree, 3 for multivariate simple tree

 SimulationTime <- 100
 blank_matrix <- matrix(1)
 tree.control=rpart.control(cp=1e-4, maxsurrogate=5, usesurrogate=2)



 taskId <- Sys.getenv("SLURM_ARRAY_TASK_ID")

 taskId <- as.numeric(taskId)

 #taskId <- 650
 
 tempId <- ceiling(taskId / 5)
 I <-  I_seq[taskId - (tempId - 1) * 5] 
 tempId2 <- ceiling(tempId  / 4)
 sigmae <- sigmae_seq[tempId - (tempId2 - 1) * 4 ]
 tempId3 <- ceiling(tempId2  / 4)
 Ti <- Ti_seq[tempId2 - (tempId3 - 1) * 4 ]
 structure <- structure_seq[tempId3]

#structure <- 3
#I <- 200
#Ti <- 10
#sigmae <- 1

 if(structure < 3){
    
    VariableNumber <- 2
    
    if(structure == 1){
       load("tree_true_1.Rdata") 
       
       FixedTree <- function(X, mu){
          Indicator <- (X < 5)
          Leaves <- cbind((Indicator[ ,1] & Indicator[ ,2]),
                          (Indicator[ ,1] & (!Indicator[ ,2])),
                          ((!Indicator[ ,1]) & Indicator[ ,3]),
                          ((!Indicator[ ,1]) & (!Indicator[ ,3])))
          return(Leaves %*% mu)
       }
       
       LeafValue <- matrix(c(10:13, 6:9), ncol=VariableNumber)
       
       p <- 7
       
       
    }else{
       load("tree_true_2.Rdata") 
       
       FixedTree <- function(X, mu){
          Indicator <- (X < 5)
          Leaves <- cbind((Indicator[ ,1] & Indicator[ ,2] & Indicator[ ,4]),
                          (Indicator[ ,1] & Indicator[ ,2] & !Indicator[ ,4]),
                          (Indicator[ ,1] & (!Indicator[ ,2]) & Indicator[ ,5]),
                          (Indicator[ ,1] & (!Indicator[ ,2]) & !Indicator[ ,5]),
                          (!Indicator[ ,1] & Indicator[ ,3]),
                          (!Indicator[ ,1] & (!Indicator[ ,3]) & Indicator[ ,6]),
                          (!Indicator[ ,1] & (!Indicator[ ,3]) & (!Indicator[ ,6])))
          return(Leaves %*% mu)
       }
       
       p <- 10
       
       LeafValue <- matrix(c(seq(6, 18, by=2), seq(4.5, 16.5, by=2)), ncol=VariableNumber)
    }
   
 }else{
    
    VariableNumber <- 5
    
    
    load("tree_true_1.Rdata") 
       
    FixedTree <- function(X, mu){
          Indicator <- (X < 5)
          Leaves <- cbind((Indicator[ ,1] & Indicator[ ,2]),
                          (Indicator[ ,1] & (!Indicator[ ,2])),
                          ((!Indicator[ ,1]) & Indicator[ ,3]),
                          ((!Indicator[ ,1]) & (!Indicator[ ,3])))
          return(Leaves %*% mu)
       }
       
    p <- 7
       
    LeafValue <- matrix(c(10:13, 9:12, 8:11, 6:9, 4:7), ncol=VariableNumber)
      
}

ResponseVariable <- paste("V", 1:VariableNumber, sep="")

        
set.seed(taskId)


results <- foreach(t=1:SimulationTime, 
                   .combine='rbind', 
                   .packages = c("nlme", "MASS","mvpart","expm", "Formula", "stats"), 
                   .errorhandling="remove") %dorng%{        
#for(t in 1:SimulationTime)  {
   N <- I * Ti
   
   # Generate predictors
   X <- matrix(runif(p * N, min=0, max=10), ncol=p)
   
   X[ ,p - 1] <- round(X[ ,p - 1])
   X[ ,p] <- 5.77 * rbinom(N, 1, 0.5)
   
   FixedValue <- FixedTree(X, LeafValue)
   
   group <- as.factor(rep(1:I, each=Ti))
   
   simb <- 0
   
   sime <- mvrnorm(N, rep(0, VariableNumber), sigmae * diag(VariableNumber)) # the i-th row of sime  =  (e_iY, e_iW)
   
   # Generate responses
   Responses <- FixedValue + sime
   
   colnames(Responses) <- ResponseVariable
   
   # Construct the training set
   
   train <- data.frame(X, group, Responses)
   
   Predictors <- paste(colnames(train)[1:p], collapse = "+")
   
   # Generate a test set with size 20I
   T_test <- 20
   N_test <- T_test * I
   
   X_test <- matrix(runif(p * N_test, min=0, max=10), ncol=p)
   
   X_test[ ,p-1] <- round(X_test[ ,p-1])
   X_test[ ,p] <- 5.77 * rbinom(N_test, 1, 0.5)
   
   FixedValue_test <- FixedTree(X_test, LeafValue)
   
   group_test <- as.factor(rep(1:I, each=T_test))
   
   e_test <- mvrnorm(N_test,rep(0, VariableNumber), sigmae * diag(VariableNumber))
   
   Responses_test <- FixedValue_test + e_test
   
   test <- data.frame(X_test, group_test, Responses_test)
   names(test) <- names(train)

   
   # Method 1 : Seperate univariate lmes

   predict_unilme_fixed <- matrix(NA, nrow=nrow(test), ncol=VariableNumber)

   predict_unilme_group <- matrix(NA, nrow=nrow(test), ncol=VariableNumber)


   for (i in 1:VariableNumber) {


      RV <- ResponseVariable[i]

      Formula_unilme <- as.formula(paste(RV, "~", Predictors, sep=""))

      .GlobalEnv$Formula_unilme  <- Formula_unilme # function lme can only look for the variables in lme.formula from the global environment.

      fit_unilme <- lme(fixed=Formula_unilme,
                        data=train,
                        random = ~ 1 | group,
                        method='REML',
                        control=lmeControl(opt = "optim"))


      predict_unilme_fixed[ ,i] <- as.numeric(predict(fit_unilme, newdata = test, level=0))
      predict_unilme_group[ ,i] <- as.numeric(predict(fit_unilme, newdata = test, level=1))


   }



   mse_unilme_fixed <- sum((predict_unilme_fixed - FixedValue_test)^2)/ nrow(test)
   mse_unilme_group <- sum((predict_unilme_group - test[ ,ResponseVariable])^2)/ nrow(test)



   # Method 2 : One multivariate lme

   value <- c(as.matrix(train[ ,ResponseVariable]))
   trait <- as.factor(rep(1:VariableNumber, each=nrow(train)))
   train_stacked <- cbind(value, train[rep(1:nrow(train), VariableNumber), ], trait)

   Formula_multilme  <- as.formula(paste("value ~ 0 + trait + ",paste("trait:", names(train)[1:p], collapse = "+")))

   .GlobalEnv$Formula_multilme  <- Formula_multilme

   fit_multilme <- lme(Formula_multilme,
                       data=train_stacked,
                       random=  ~ 0 + trait | group,
                       method='REML',
                       control=lmeControl(opt = "optim"))

   value_test <- c(as.matrix(test[ ,ResponseVariable]))
   trait_test <- as.factor(rep(1:VariableNumber, each=nrow(test)))
   test_stacked <- data.frame(value=value_test, test[rep(1:nrow(test), VariableNumber), ], trait=trait_test)

   predict_multilme_fixed <- as.numeric(predict(fit_multilme, newdata = test_stacked, level=0))
   predict_multilme_group <- as.numeric(predict(fit_multilme, newdata = test_stacked, level=1))

   mse_multilme_fixed <- sum((predict_multilme_fixed - c(FixedValue_test))^2)/ nrow(test)
   mse_multilme_group <- sum((predict_multilme_group - value_test)^2)/ nrow(test)
   
   
   
   
   
   # Method 3 : Single multivariate regression tree
   
   formula <- as.formula(paste(paste(ResponseVariable, collapse = "+"), "~", Predictors))
   
   fit_multitree <- multitree(formula, data=train)
   
   predict_multitree <- predict.multitree(fit_multitree, newdata = test)
   
   emse_multitree <- sum((predict_multitree - FixedValue_test)^2)/ nrow(test)
   mse_multitree <- sum((predict_multitree - test[ ,ResponseVariable])^2)/ nrow(test)
   
   recover_multitree <- is.sametree.rpart(tree_true, fit_multitree$Tree)
   
   #tree_true <- fit_multitree
   
   #save(tree_true, file="tree_true_2.Rdata")
   
   
   # Method 4: Multiple univariate REEM-tree 
   
   predict_uniREEM_fixed <- matrix(NA, nrow=nrow(test), ncol=VariableNumber)
   
   predict_uniREEM_group <- matrix(NA, nrow=nrow(test), ncol=VariableNumber)
   
   recover_uniREEM <- rep(NA, VariableNumber)
   
   RandomEffects_uniREEM <-matrix(NA, nrow=I, ncol=VariableNumber)
   
   for (i in 1:VariableNumber) {
      
      RV <- ResponseVariable[i]
      
      fit.REEM <- REEMtree(as.formula(paste(RV, "~", Predictors)), 
                         data=train, 
                         random= ~ 1 | group, 
                         subset=NULL, 
                         initialRandomEffects=0,
                         ErrorTolerance=0.001, 
                         MaxIterations=1000, 
                         verbose=FALSE, 
                         tree.control=rpart.control(),
                         cv=TRUE, 
                         cpmin = 0.01, 
                         no.SE = 1,
                         lme.control=lmeControl(opt = "optim"), 
                         method="REML", 
                         correlation=NULL
                       ) 
      
      tree_uniREEM <- fit.REEM$Tree
      
      # whether uniREEM separately recovers the true tree struture
      recover_uniREEM[i] <- is.sametree.rpart(tree_true, tree_uniREEM)
      
      TreePrediction <- predict(tree_uniREEM, test)

      
      # Prediciton
      predict_uniREEM_fixed[ ,i] <- TreePrediction
      predict_uniREEM_group[ ,i] <- TreePrediction  +  fit.REEM$RandomEffects[test$group, ]

      RandomEffects_uniREEM[ ,i] <- as.matrix(fit.REEM$RandomEffects)
      
   }
   

   mse_uniREEM_fixed <- sum((predict_uniREEM_fixed - FixedValue_test)^2)/ nrow(test)
   mse_uniREEM_group <- sum((predict_uniREEM_group - test[ ,ResponseVariable])^2)/ nrow(test)
   
   RandomEffects_EMSE_uniREEM <- mean((as.matrix(RandomEffects_uniREEM - simb))^2)
   
  
   
   # 5. Proposed Method: multiREEMtree_1se_cov  

   
   formula <- as.formula(paste(paste(ResponseVariable, collapse = "+"), "~", Predictors))
   
   
   fit_multiREEM_1se_cov <- multiREEMtree(formula,
                                          data=train,
                                          random= ~ 1 | group,
                                          subset=NULL,
                                          stand.y = "covariance",
                                          tree.xv='1se',
                                          tree.xval=10,
                                          tree.control=rpart.control(cp=1e-4, maxsurrogate=5, usesurrogate=2),
                                          lme.algorithm='nlme',
                                          correlation=NULL)
   
   
   
   tree_multiREEM_1se_cov <- fit_multiREEM_1se_cov$Tree
   
   recover_multiREEMtree_1se_cov <- is.sametree.rpart(tree_true, tree_multiREEM_1se_cov)
   
   predict_multiREEMtree_1se_cov_fixed  <- predict.multiREEMtree(fit_multiREEM_1se_cov, newdata=test, level=0)
   predict_multiREEMtree_1se_cov_group  <- predict.multiREEMtree(fit_multiREEM_1se_cov, newdata=test, level=1)
   
   mse_multiREEMtree_1se_cov_fixed <- sum((predict_multiREEMtree_1se_cov_fixed - FixedValue_test)^2) / nrow(test)
   mse_multiREEMtree_1se_cov_group <- sum((predict_multiREEMtree_1se_cov_group - test[ ,ResponseVariable])^2) / nrow(test)
   
   RandomEffects_EMSE_multiREEMtree_1se_cov <- mean((as.matrix(fit_multiREEM_1se_cov$RandomEffects$group - simb))^2)
   

   
   # 6. Proposed Method: multiREEMtree_min_cov  
   

     
   fit_multiREEM_min_cov <-  multiREEMtree(formula,
                                           data=train,
                                           random= ~ 1 | group,
                                           subset=NULL,
                                           stand.y = "covariance",
                                           tree.xv='min',
                                           tree.xval=10,
                                             tree.control=rpart.control(cp=1e-4, maxsurrogate=5, usesurrogate=2),
                                           lme.algorithm='nlme',
                                           correlation=NULL)
     
     
   tree_multiREEM_min_cov <- fit_multiREEM_min_cov$Tree
   
   recover_multiREEMtree_min_cov <- is.sametree.rpart(tree_true, tree_multiREEM_min_cov)
   
   predict_biREEMtree_min_cov_fixed  <- predict.multiREEMtree(fit_multiREEM_min_cov, newdata=test, level=0)
   predict_biREEMtree_min_cov_group  <- predict.multiREEMtree(fit_multiREEM_min_cov, newdata=test, level=1)
   
   
   mse_multiREEMtree_min_cov_fixed <- sum((predict_biREEMtree_min_cov_fixed - FixedValue_test)^2) / nrow(test)
   mse_multiREEMtree_min_cov_group <- sum((predict_biREEMtree_min_cov_group - test[ ,ResponseVariable])^2) / nrow(test)
   
   RandomEffects_EMSE_multiREEMtree_min_cov <- mean((as.matrix(fit_multiREEM_min_cov$RandomEffects$group - simb))^2)
   
   
   # 7. Proposed Method: multiREEMtree_1se_marg
   
   fit_multiREEM_1se_marg  <- multiREEMtree(formula,
                                            data=train,
                                            random= ~ 1 | group,
                                            subset=NULL,
                                            stand.y = "marginal",
                                            tree.xv='1se',
                                            tree.xval=10,
                                            tree.control=rpart.control(cp=1e-4, maxsurrogate=5, usesurrogate=2),
                                            lme.algorithm='nlme',
                                            correlation=NULL)
   
   
   tree_multiREEM_1se_marg <- fit_multiREEM_1se_marg$Tree
   
   recover_multiREEMtree_1se_marg <- is.sametree.rpart(tree_true, tree_multiREEM_1se_marg)
   
   predict_multiREEMtree_1se_marg_fixed  <- predict.multiREEMtree(fit_multiREEM_1se_marg, newdata=test, level=0)
   predict_multiREEMtree_1se_marg_group  <- predict.multiREEMtree(fit_multiREEM_1se_marg, newdata=test, level=1)
   
   
   mse_multiREEMtree_1se_marg_fixed <- sum((predict_multiREEMtree_1se_marg_fixed - FixedValue_test)^2) / nrow(test)
   mse_multiREEMtree_1se_marg_group <- sum((predict_multiREEMtree_1se_marg_group - test[ ,ResponseVariable])^2) / nrow(test)
   

   RandomEffects_EMSE_multiREEMtree_1se_marg <- mean((as.matrix(fit_multiREEM_1se_marg$RandomEffects$group - simb)^2))
   

   
   # 8. Proposed Method: multiREEMtree_min_marg
   
   fit_multiREEM_min_marg  <- multiREEMtree(formula,
                                            data=train,
                                            random= ~ 1 | group,
                                            subset=NULL,
                                            stand.y = "marginal",
                                            tree.xv='min',
                                            tree.xval=10,
                                            tree.control=rpart.control(cp=1e-4, maxsurrogate=5, usesurrogate=2),
                                            lme.algorithm='nlme',
                                            correlation=NULL)
   
   tree_multiREEM_min_marg <- fit_multiREEM_min_marg$Tree
   
   recover_multiREEMtree_min_marg <- is.sametree.rpart(tree_true, tree_multiREEM_min_marg)
   
   predict_multiREEMtree_min_marg_fixed  <- predict.multiREEMtree(fit_multiREEM_min_marg, newdata=test, level=0)
   predict_multiREEMtree_min_marg_group  <- predict.multiREEMtree(fit_multiREEM_min_marg, newdata=test, level=1)
   
   
   mse_multiREEMtree_min_marg_fixed <- sum((predict_multiREEMtree_min_marg_fixed - FixedValue_test)^2) / nrow(test)
   mse_multiREEMtree_min_marg_group <- sum((predict_multiREEMtree_min_marg_group - test[ ,ResponseVariable])^2) / nrow(test)
   
   RandomEffects_EMSE_multiREEMtree_min_marg <- mean((as.matrix(fit_multiREEM_min_marg$RandomEffects$group - simb)^2))
   
   #if(t %% 100 == 0){
   #   save(blank_matrix, file=paste("eyes//eyes, I=", I, ", T=", Ti, ", cor=", sigmab, ", noise=", sigmae, ", run=", t, ".RData", sep=""))
   #}
   
   
   
   # results <- rbind(results, c(mse_unilme_fixed, mse_multilme_fixed, mse_multitree, mse_uniREEM_fixed,
   #                             mse_multiREEMtree_1se_fixed, mse_multiREEMtree_min_cov_fixed, mse_multiREEMtree_1se_marg_fixed,
   #                             mse_unilme_group, mse_multilme_group, mse_uniREEM_group, 
   #                             mse_multiREEMtree_1se_group, mse_multiREEMtree_min_cov_group, mse_multiREEMtree_1se_marg_group,
   #                             recover_multitree, recover_uniREEM_Y, recover_uniREEM_W, 
   #                             recover_multiREEMtree_1se, recover_multiREEMtree_min_cov, recover_multiREEMtree_1se_marg))
   
   return(c(mse_unilme_fixed, mse_multilme_fixed, emse_multitree, mse_uniREEM_fixed,
            mse_multiREEMtree_1se_cov_fixed, mse_multiREEMtree_min_cov_fixed, mse_multiREEMtree_1se_marg_fixed, mse_multiREEMtree_min_marg_fixed,
            mse_unilme_group, mse_multilme_group, mse_multitree, mse_uniREEM_group, 
            mse_multiREEMtree_1se_cov_group, mse_multiREEMtree_min_cov_group, mse_multiREEMtree_1se_marg_group, mse_multiREEMtree_min_marg_group,
            recover_multitree, recover_uniREEM, 
            recover_multiREEMtree_1se_cov, recover_multiREEMtree_min_cov, recover_multiREEMtree_1se_marg, recover_multiREEMtree_min_marg,
            RandomEffects_EMSE_uniREEM, RandomEffects_EMSE_multiREEMtree_1se_cov,  RandomEffects_EMSE_multiREEMtree_min_cov,
            RandomEffects_EMSE_multiREEMtree_1se_marg, RandomEffects_EMSE_multiREEMtree_min_marg))
            

            
}

        
colnames(results) <- c("mse_unilme_fix", "mse_multilme_fix", "emse_multitree", "mse_uniREEM_fixed",
                               "mse_1se_cov_fixed", "mse_min_cov_fixed", "mse_1se_marg_fixed","mse_min_marg_fixed",
                               "mse_unilme_group", "mse_multilme_group", "mse_multitree", "mse_uniREEM_group", 
                               "mse_1se_cov_group", "mse_min_cov_group", "mse_1se_marg_group", "mse_min_marg_group",
                               "recover_multitree", paste("recover_uniREEM", ResponseVariable, sep="_"), 
                               "recover_1se_cov", "recover_min_cov", "recover_1se_marg", "recover_min_marg",
                       "RE_uniREEM", "RE_1se_cov", "RE_min_cov", "RE_1se_marg", "RE_min_marg")
save(results, file=paste("results_no_random//results, I=", I, ", T=", Ti, ", noise=", sigmae, ", structure=", structure, ".RData", sep=""))
        
#      }
 #   }
#  }
#} 

stopCluster(cl)




