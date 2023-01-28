#' Wrapper function for calling of mvpart function
#'
#' Wrapper function for calling of mvpart function
#'
#' @param formula A formula of the form \code{y1 + ... + ym ~ x1 + ... + xn}, where \code{y1,...,ym} are response
#' variables and \code{x1,...,xn} are predictors.
#' @param data A data frame that contains all of the variables named in the formula.
#' @param subset  A logical vector showing the subset of data to use for fitting. Defaults to be a vector of TRUEs.
#' @param stand.y An option to standardize the response variables, should be one of "marginal" (standardize each response variable by subtracting
#' its mean and dividing by its sample standard deviation), "covariance"(standardize the set of response variables using the vector of sample
#' means and sample covariance matrix, i.e., convert y to y * (cov(y)^(-1/2)), or "none" (no standardization). Defaults to be "marginal".
#' @param tree.xv Selection of tree by cross-validation: \code{"1se"} (default) - gives the simplest tree that has cross-validation error within
#' one SE of the minimum value using k-fold cross-validation, with k given by tree.xval; \code{"min"} - the tree with minimum error
#' of k-fold cross-validation, with k given by tree.xval; \code{"none"} - no cross-validation.
#' @param tree.xval Number of cross-validation folds or vector defining cross-validation groups. Defaults to be 10.
#' Could also be \code{"LOOCV"}, which will use leave-one-out cross-validation.
#' @param tree.control A list of control values used to control the \link[mvpart]{rpart}; See \link[mvpart]{rpart.control}
#' for details. Defaults to be \code{rpart.control(cp=0.01, maxsurrogate=5, usesurrogate=2)}.
#'
#' @return A multitree object, which includes.
#' \itemize{
#'  \item Tree - The fitted multivariate tree.
#'  \item ResponseVariables - The response variables \code{y1,...,ym} defined by \code{formula}.
#'  \item Predictors - The predictors \code{x1,...,xn} defined by \code{formula}.
#'  \item data - The data set.
#'  \item Formula - The parameter \code{formula} in the function.
#'  \item Subset - The parameter \code{Subset} in the function.
#'  \item tree.control- \code{tree.control} in the function.
#'  \item xv - The parameter \code{tree.xv} in the function.
#'  \item xval - The parameter \code{tree.xval} in the function.
#'  \item Y_original - The original values of the response variables.
#'  \item stand.y - The parameter \code{stand.y} in the function.
#' }
#'
#' @import mvpart
#' @import Formula
#' @importFrom stats as.formula terms cov sd
#'
#' @export

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
  
  result <- list(Tree = tree, Formula = formula, ResponseVariables = ResponseVariables, 
                 Predictors = Predictors, data = data, Subset = SubsetVector, 
                 xv = tree.xv, xval = tree.xval, tree.control = tree.control, Y_original = OriginalResponse, 
                 stand.y = stand.y)
  class(result) <- "multitree"
  return(result)
}

