#' @title Selection of best LDA model(s) for use in time series
#'
#' @param lda_models LDA model output
#' @param LDA_eval function name for evaluation of the LDA models
#' @param LDA_selector function name for selecting the LDA model(s) 
#' @param ... additional arguments to be passed to subfunctions
#' @return the output of the selected LDA model(s)
#'
#' @export
#'
LDA_select <- function(lda_models = NULL, LDA_eval = quote(AIC),  
                       LDA_selector = quote(min), ...){

  lda_eval <- sapply(lda_models, LDA_eval) %>%
              matrix(ncol = 1)
  lda_selected <- apply(lda_eval, 2, LDA_selector) 
  which_selected <- which(lda_eval %in% lda_selected)
  selected_lda <- lda_models[[which_selected]]
  out <- selected_lda
  return(out)
}

#' @title Prep the data inputs to the MTS function
#'
#' @param lda_models selected LDA model(s)
#' @param document_covariate_matrix matrix of documents (rows) by covariates
#'   (columns)
#' @return data matrix ready for MTS analyses
#'
#' @export
#'
MTS_prep <- function(lda_models = NULL, document_covariate_matrix = NULL){

  nmods <- length(lda_models)
  if(("LDA_list" %in% class(lda_models)) == FALSE){
    lda_models <- list(lda_models)
    nmods <- 1
  }
  out <- vector("list", length = nmods)
  for(i in 1:nmods){
    out[[i]] <- data.frame(document_covariate_matrix)
    out[[i]]$gamma <- lda_models[[i]]@gamma
  }
  return(out)
}

#' @title Conduct a set of MTS analyses
#'
#' @param prepped_data data prepped for each of the MTS analyses
#' @param formula vector of formulas for the continuous change
#' @param nchangepoints vector of the number of change points to include in 
#'   the model
#' @param vector of weights for each document
#' @param ... additional arguments to be passed to subfunctions
#' @return 
#'
#' @export
#'
MTS_set <- function(data = NULL, formula, nchangepoints,  
                    weights = rep(1, nrow(data[[1]])), ...){

  ldas <- 1:length(data)
  if(is(formula, "formula")){
    formula <- c(formula)
  }
  mods <- expand.grid(lda = ldas, formula = formula, 
                      nchangepoints = nchangepoints, stringsAsFactors = FALSE)
  nmods <- nrow(mods)
  out <- vector("list", nmods)
  for(i in 1:nmods){
    out[[i]] <- LDATS::MTS(data[[mods$lda[i]]], mods$formula[i][[1]], 
                           mods$nchangepoints[i], weights, ...)
  }
  return(out)
}