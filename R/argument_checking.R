

#' @title Check that a control list is proper
#' 
#' @description Check that a list of controls is of the right class.
#'   
#' @param control Control list to evaluate.
#' 
#' @param eclass Expected class of the list to be evaluated.
#'
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  check_control(list())
#'
#' @export
#'
check_control <- function(control, eclass = "list"){
  if (!(eclass %in% class(control))){
    stop(paste0("control is not a ", eclass))
  }
  return()
}


#' @title Check that document term table is proper
#' 
#' @description Check that the table of observations is conformable to
#'   a matrix of integers.
#'   
#' @param document_term_table Table of observation count data (rows: 
#'   documents, columns: terms. May be a class \code{matrix} or 
#'   \code{data.frame} but must be conformable to a matrix of integers,
#'   as verified by \code{\link{check_document_term_table}}. 
#' 
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  data(rodents)
#'  check_document_term_table(rodents$document_term_table)
#'
#' @export
#'
check_document_term_table <- function(document_term_table){
  document_term_table_m <- as.matrix(document_term_table)
  if(any(document_term_table_m %% 1 != 0)){
    dtt <- "document_term_table"
    msg <- paste0(dtt, " is not conformable to a matrix of integers")
    stop(msg)
  }
  return()
}

#' @title Check that topics vector is proper
#' 
#' @description Check that the vector of numbers of topics is conformable to
#'   integers greater than 1.
#'   
#' @param topics Vector of the number of topics to evaluate for each model.
#'   Must be conformable to \code{integer} values.
#'
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  check_topics(2)
#'
#' @export
#'
check_topics <- function(topics){
  if (!is.numeric(topics) || any(topics %% 1 != 0)){
    stop("topics vector must be integers")
  }
  if (any(topics < 2)){
    stop("minimum number of topics currently allowed is 2")
  }
  return()
}

#' @title Check that nseeds value or vector is proper
#' 
#' @description Check that the vector of numbers of seeds is conformable to
#'   integers greater than 0.
#'   
#' @param nseeds \code{integer} number of seeds (replicate starts) to use for 
#'   each value of \code{topics} in the LDAs. Must be conformable to a  
#'   positive \code{integer} value.
#' 
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  check_nseeds(1)
#'  check_nseeds(2)
#'
#' @export
#'
check_nreps <- function(nreps){
  if (!is.numeric(nreps) || any(nreps %% 1 != 0)){
    stop("nreps must be integer-conformable")
  }
  return()
}

check_nchangepoints <- function(nchangepoints){
  if (!is.numeric(nchangepoints) || any(nchangepoints %% 1 != 0)){
    stop("nchangepoints must be integer-valued")
  }
  if (any(nchangepoints < 0)){
    stop("nchangepoints must be non-negative")
  }
  return()
}

check_weights <- function(weights){
  if(is.logical(weights)){
    if(weights){
      return()
    } else{
      stop("if logical, weights need to be TRUE")
    }   
  }
  if(!is.null(weights)){
    if (!is.numeric(weights)){
      stop("weights vector must be numeric")
    }
    if (any(weights <= 0)){
      stop("weights must be positive")
    }
    if (round(mean(weights)) != 1){
      warning("weights should have a mean of 1, fit may be unstable")
    }
  }
  return()
}


check_LDA_models <- function(LDA_models){
  if(("LDA_set" %in% class(LDA_models)) == FALSE){
    if(is(LDA_models, "LDA") == FALSE){
      stop("LDA_models is not an LDA object or LDA_set object")
    }
  }
  return()
}

check_document_covariate_table <- function(document_covariate_table, 
                                           LDA_models = NULL,
                                           document_term_table = NULL){
  dct_df <- tryCatch(data.frame(document_covariate_table),
                     warning = function(x){NA}, error = function(x){NA})
  if(is(LDA_models, "LDA")){
    LDA_models <- c(LDA_models)
    class(LDA_models) <- c("LDA_set", "list")
  }
  if (length(dct_df) == 1 && is.na(dct_df)){
    stop("document_covariate_table is not conformable to a data frame")
  }
  if (!is.null(LDA_models)){
    if (nrow(data.frame(document_covariate_table)) != 
        nrow(LDA_models[[1]]@gamma)){
      stop("number of documents in covariate table is not equal to number of 
        documents observed")
    }
  } else if (!is.null(document_term_table)){
    if (nrow(data.frame(document_covariate_table)) != 
        nrow(data.frame(document_term_table))){
      stop("number of documents in covariate table is not equal to number of 
        documents observed")
    }
  }
  return()
}

check_timename <- function(document_covariate_table, timename){
  if (!("character" %in% class(timename))){
    stop("timename is not a character value")
  }
  if (length(timename) > 1){
    stop("timename can only be one value")
  }
  covariate_names <- colnames(document_covariate_table)
  if ((timename %in% covariate_names) == FALSE){
    stop("timename not present in document covariate table")
  }
  time_covariate <- document_covariate_table[ , timename]
  if (!(is.Date(time_covariate)) & 
      (!is.numeric(time_covariate) || !all(time_covariate %% 1 == 0))){
    stop("covariate indicated by timename is not an integer or a date")
  }
  return()
}


check_formulas <- function(formulas, document_covariate_table, 
                           control = list()){
  check_document_covariate_table(document_covariate_table)
  check_control(control)
  control <- do.call("TS_control", control)
  # response <- control$response
  dct <- document_covariate_table
  if (!is(formulas, "list")) {
    if (is(formulas, "formula")) {
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain formula(s)")
    }
  } else if (!all(vapply(formulas, is, TRUE, "formula"))) {
      stop("formulas does not contain all formula(s)")
  }
  resp <- unlist(lapply(lapply(formulas, terms), attr, "response"))
  pred <- unlist(lapply(lapply(formulas, terms), attr, "term.labels"))
  if (any(resp != 0)) {
    stop("formula inputs should not include response variable")
  }
  if (!all(pred %in% colnames(dct))) {
    misses <- pred[which(pred %in% colnames(dct) == FALSE)]
    mis <- paste(misses, collapse = ", ")
    stop(paste0("formulas include predictors not present in data: ", mis))
  }
  return()
}




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

