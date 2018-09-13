#' @title Conduct a set of Time Series analyses on LDA models
#'
#' @description This is a wrapper function that expands the main Time Series
#'   analyses function (\code{MTS}) across the LDA models and the Time Series 
#'   models, with respect to both continuous time formulae and the number
#'   of discrete changepoints. This function allows direct passage of the
#'   control parameters for the parallel tempering MCMC through to the 
#'   main Time Series function, \code{MTS}, via the \code{ptMCMC_controls}
#'   argument  
#'
#' @param lda_models List of LDA models (class \code{LDA_list}) or a singular
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
#' @param timename Character name of the column in the 
#'   \code{document_covariate_table} that contains the time index to use
#'   for assignment of the changepoints (corresponding to the vector 
#'   \strong{\eqn{t}} in the mathematical description of the model). 
#'
#' @param formula Vector of \code{formula}(s) for the continuous change. Any 
#'   predictor variable included in a formula must also be a column in the
#'   \code{document_covariate_table}. Each element (formula) in the vector
#'   is evaluated for each number of change points and each LDA model.
#'
#' @param changepoints Vector of integers corresponding to the number of 
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
#' @param ptMCMC_controls Class \code{ptMCMC_controls} list, holding control 
#'   parameters for the parallel tempering Markov Chain Monte Carlo (ptMCMC)
#'   fitting of the MTS models.
#'
#' @return the set of results from MTS
#'
#' @export
#'
TS_set_on_LDA <- function(lda_models, document_covariate_table, timename,
                          formula = ~ 1, changepoints = 0, weights = NULL, 
                          ptMCMC_controls = ptMCMC_controls_list()){

  lda_models <- check_LDA_models(lda_models)
  check_document_covariate_table(document_covariate_table, lda_models)
  check_timename(document_covariate_table, timename)
  formula <- check_formula(formula, document_covariate_table)  
  check_changepoints(changepoints)
  check_weights(weights)

  mods <- expand_TS(lda_models, formula, changepoints)
  nmods <- nrow(mods)
  out <- vector("list", nmods)
  for(i in 1:nmods){

  }
  return(out)
}

check_changepoints <- function(changepoints){
# to do: verify that the input is a vector of integers

}

check_weights <- function(weights){
# to do: verify that the input is a vector of numeric values, should be [0,1]
}

expand_TS <- function(lda_models, formula, changepoints){
  nmods <- length(lda_models)
  expand.grid(lda = 1:nmods, formula = formula, 
                      nchangepoints = changepoints, stringsAsFactors = FALSE)
}

check_formula <- function(formula, document_covariate_table){
  if(!is(formula, "vector")){
    if(is(formula, "formula")){
      formula <- c(formula)
    } else{
      stop("formula is not a formula")
    }
  } else{ 
    if (!all(unlist(lapply(formula, is, "formula")))){
      stop("formula is not a vector of formulas")
    }
  }

  # to do: verify that all predictors in the formulas are in the table

  formula
}


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

check_document_covariate_table <- function(document_covariate_table, 
                                           lda_models){
  dct_df <- tryCatch(data.frame(document_covariate_table),
                     warning = function(x){NA}, error = function(x){NA})
  if (length(dct_df) == 1 && is.na(dct_df)){
    stop("document_covariate_table is not conformable to a data frame")
  }
  if (nrow(document_covariate_matrix) != nrow(lda_models[[1]]@gamma)){
    stop("number of documents in covariate table is not equal to number of 
      documents observed")
  }
}


check_lda_models <- function(lda_models){
  if(("LDA_list" %in% class(lda_models)) == FALSE){
    if(is(lda_models, "LDA") == TRUE){
      lda_models <- list(lda_models)
      class(lda_models) <- c("LDA_list", "list")
    } else{
      stop("lda_models is not an LDA object or LDA_list object")
    }
  }
  lda_models
}


ptMCMC_controls_list <- function(){
  out <- list()
  class(out) <- c("ptMCMC_controls", "list")
  out
}