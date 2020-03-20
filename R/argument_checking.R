#' @title Check that a formula is proper
#' 
#' @description Check that \code{formula} is actually a 
#'   \code{\link[stats]{formula}} and that the
#'   response and predictor variables are all included in \code{data}.
#'   
#' @param formula \code{formula} to evaluate.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. 
#' 
#' @return An error message is thrown if \code{formula} is not proper,
#'   else \code{NULL}.
#' 
#' @export
#'
check_formula <- function(data, formula){

  if (!is(formula, "formula")){
    stop("formula does not contain a single formula")
  }

  respLoc <- attr(terms(formula), "response")
  if (respLoc == 0){
    stop("formula inputs should include response variable")
  }  

  resp <- as.character(attr(terms(formula), "variables"))[-1][respLoc]
  pred <- attr(terms(formula), "term.labels")
  if (!resp %in% colnames(data)){
    stop("formula includes response not present in data")
  }
  if (!all(pred %in% colnames(data))){
    misses <- pred[which(pred %in% colnames(data) == FALSE)]
    mis <- paste(misses, collapse = ", ")
    stop(paste0("formula includes predictors not present in data: ", mis))
  }
  return()
}

