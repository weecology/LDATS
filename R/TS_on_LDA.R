#' @title Conduct a set of Time Series analyses on a set of LDA models
#'
#' @description This is a wrapper function that expands the main Time Series
#'   analyses function (\code{\link{TS}}) across the LDA models and the Time
#'   Series models, with respect to both continuous time formulas and the 
#'   number of discrete changepoints. This function allows direct passage of
#'   the control parameters for the parallel tempering MCMC through to the 
#'   main Time Series function, \code{\link{TS}}, via the 
#'   \code{ptMCMC_controls} argument.
#'
#' @param LDA_models List of LDA models (class \code{LDA_list}) or a singular
#'   LDA model (class \code{LDA}).
#'
#' @param document_covariate_table Document covariate table (rows:
#'   documents (\code{M}), columns: time index and covariate options). 
#'   Every model needs a covariate to describe the time value for each
#'   document (in whatever relevant units), whose name in the table is input
#'   via \code{timename}, that dictates the application of the changepoints. 
#'   In addition, the table needs to include as columns all covariates named 
#'   within the specific models described via the argument \code{formula} 
#'   (if desired). Must be a conformable to a data table. 
#'
#' @param formulas Vector of \code{formula}(s) for the continuous change. Any 
#'   predictor variable included in a formula must also be a column in the
#'   \code{document_covariate_table}. Each element (formula) in the vector
#'   is evaluated for each number of change points and each LDA model.
#'
#' @param nchangepoints Vector of integers corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. Each element 
#'   (number of change points) in the vector is used to dictate the
#'   segementation of the data  for each continuous model and each LDA model.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}.
#'
#' @return The set of results from TS.
#'
#' @export
#'
TS_on_LDA <- function(LDA_models, document_covariate_table, formulas = ~ 1, 
                      nchangepoints = 0, weights = NULL, 
                      control = TS_controls_list()){

  LDA_models <- check_LDA_models(LDA_models)
  check_document_covariate_table(document_covariate_table, LDA_models)
  check_timename(document_covariate_table, control$timename)
  formulas <- check_formulas(formulas, document_covariate_table, control)  
  check_nchangepoints(nchangepoints)
  check_weights(weights)

  mods <- expand_TS(LDA_models, formulas, nchangepoints)
  nmods <- nrow(mods)
  out <- vector("list", nmods)
  for(i in 1:nmods){
    formula_i <- mods$formula[[i]]
    nchangepoints_i <- mods$nchangepoints[i]
    data_i <- document_covariate_table
    data_i$gamma <- LDA_models[[mods$LDA[i]]]@gamma
    out[[i]] <- TS(data_i, formula_i, nchangepoints_i, weights, control)
  }
  return(out)
}

#' @title Expand the TS models needed across the factorial combination of
#'   LDA models, continuous formulas, and number of change points
#' 
#' @description Expand the completely crossed combination of model inputs: 
#'   LDA model results, continuous formulas
#'   
#' @param LDA_models \code{LDA_list}-class object of LDA models.
#' 
#' @param formulas Vector of the continuous formulas. 
#'
#' @param nchangepoints Vector of the number of changepoints.
#'
#' @return Expanded table of the three values: [1] the LDA model (indicated
#'   as a numeric element reference to the \code{LDA_list} object), [2] the 
#'   continuous formula, and [3] the number of changepoints.
#' 
#' @export
#'
expand_TS <- function(LDA_models, formulas, nchangepoints){
  nmods <- length(LDA_models)
  mods <- 1:nmods
  out <- expand.grid(mods, formulas, nchangepoints, stringsAsFactors = FALSE)
  colnames(out) <- c("LDA", "formula", "nchangepoints") 
  out
}

#' @title Verify that nchangepoints vector is proper
#' 
#' @description Verify that the vector of numbers of changepoints is 
#'   conformable to integers greater than 1.
#'   
#' @param nchangepoints Vector of the number of changepoints to evaluate.
#'
#' @return Nothing.
#' 
#' @export
#'
check_nchangepoints <- function(nchangepoints){
  if (!is.numeric(nchangepoints) || any(nchangepoints %% 1 != 0)){
    stop("nchangepoints must be integer-valued")
  }
}

#' @title Verify that weights vector is proper
#' 
#' @description Verify that the vector of document weights is numeric
#'   and inform the user if weights are outside the optimal range: 
#'   \eqn{(0,1]}.
#'   
#' @param weights Vector of the document weights to evaluate.
#'
#' @return Nothing.
#' 
#' @export
#'
check_weights <- function(weights){
  if(!is.null(weights)){
    if (!is.numeric(weights)){
      stop("weights vector must be numeric")
    }
    if (min(weights) <= 0 | max(weights) > 1){
      ideal <- "weights should be scaled to (0,1]; "
      wrange <- paste0("min: ", min(weights), ", max: ", max(weights), "; ")
      warn <- "fit may be unstable"
      warning(paste0(ideal, wrange, warn))
    }
  }
}

#' @title Verify that LDA model input is proper
#' 
#' @description Verify that the \code{LDA_models} input is, in fact, LDA 
#'   models or a singular LDA model. If there is only one model, convert it
#'   from a singular class \code{LDA} list to a length-1 class \code{LDA_list}
#'   class and return it.
#'   
#' @param LDA_models List of LDA models or singular LDA model to evaluate.
#'
#' @return \code{LDA_models} as an \code{LDA_list}-class object.
#' 
#' @export
#'
check_LDA_models <- function(LDA_models){
  if(("LDA_list" %in% class(LDA_models)) == FALSE){
    if(is(LDA_models, "LDA") == TRUE){
      LDA_models <- list(LDA_models)
      class(LDA_models) <- c("LDA_list", "list")
    } else{
      stop("LDA_models is not an LDA object or LDA_list object")
    }
  }
  LDA_models
}

#' @title Verify that the document covariate table is proper
#' 
#' @description Verify that the table of document-level covariates is 
#'   conformable to a data frame and of the right size (correct number of 
#'   documents) for the document-topic output from the LDA models.
#'   
#' @param document_covariate_table Document covariate table to evaluate.
#'
#' @param LDA_models Reference LDA model list (class \code{LDA_list}) that 
#'   includes as its first element a properly fitted \code{LDA} model with 
#'   a \code{gamma} slot with the document-topic distribution. 
#'
#' @return Nothing.
#' 
#' @export
#'
check_document_covariate_table <- function(document_covariate_table, 
                                           LDA_models){
  dct_df <- tryCatch(data.frame(document_covariate_table),
                     warning = function(x){NA}, error = function(x){NA})
  if (length(dct_df) == 1 && is.na(dct_df)){
    stop("document_covariate_table is not conformable to a data frame")
  }
  if (nrow(document_covariate_table) != nrow(LDA_models[[1]]@gamma)){
    stop("number of documents in covariate table is not equal to number of 
      documents observed")
  }
}

#' @title Verify that the time vector is proper
#' 
#' @description Verify that the vector of time values is included in the 
#'   document covariate table and that it is either numeric or a date.
#'   
#' @param document_covariate_table Document covariate table used to query
#'   for the time column.
#'
#' @param timename Column name for the time variable to evaluate.
#'
#' @return Nothing.
#' 
#' @export
#'
check_timename <- function(document_covariate_table, timename){
  covariate_names <- colnames(document_covariate_table)
  if ((timename %in% covariate_names) == FALSE){
    stop("timename not present in document covariate table")
  }
  time_covariate <- document_covariate_table[ , timename]
  if (!(is.numeric(time_covariate)) & !(is.Date(time_covariate))){
    stop("covariate indicated by timename is not numeric or temporal")
  }
}

#' @title Verify that formulas vector is proper and append the response 
#'   variable
#' 
#' @description Verify that the vector of formulas is actually formatted
#'   as a vector formula objects and that the predictor variables are all 
#'   included in the document covariate table. Add the response variable
#'   to each of the formulas.
#'   
#' @param formulas Vector of the formulas to evaluate.
#'
#' @param document_covariate_table Document covariate table used to evaluate
#'   the availability of the data required by the formula inputs.
#'
#' @param control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model.
#'
#' @return Updated \code{formulas}.
#' 
#' @export
#'
check_formulas <- function(formulas, document_covariate_table, control){
  response <- control$response
  dct <- document_covariate_table
  if(!is(formulas, "vector")){
    if(is(formulas, "formula")){
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain a formula")
    }
  } else{ 
    if (!all(unlist(lapply(formulas, is, "formula")))){
      stop("formulas is not a vector of formulas")
    }
  }
  resp <- unlist(lapply(lapply(formulas, terms), attr, "response"))
  pred <- unlist(lapply(lapply(formulas, terms), attr, "term.labels"))
  if (any(resp != 0)){
    stop("formula inputs should not include response variable")
  }
  if (!all(pred %in% colnames(dct))){
    misses <- pred[which(pred %in% colnames(dct) == FALSE)]
    mis <- paste(misses, collapse = ", ")
    stop(paste0("formulas include predictors not present in data: ", mis))
  }
  out <- formulas
  for (i in 1:length(formulas)){
    tformula <- paste(as.character(formulas[[i]]), collapse = "")
    out[[i]] <- as.formula(paste(response, tformula))
  }
  out
}