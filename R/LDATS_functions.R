#' @title Two-stage LDA-time series analysis
#'
#' @param document_term_matrix matrix of documents (rows) by terms (columns)
#' @param document_covariate_matrix matrix of documents (rows) by covariates
#'   (columns)
#' @param formula vector of formulas for the continuous change
#' @param nchangepoints vector of the number of change points to include in 
#'   the model
#' @param ... additional arguments to be passed to subfunctions
#' @return (currently) the prepped data for the TS model
#'
#' @export
#'
LDA_TS <- function(document_term_matrix = NULL, 
                   document_covariate_matrix = NULL, 
                   formula = "1", nchangepoints = 1, ...){

  wts <- doc_weights(document_term_matrix)

  out <- LDATS::LDA(data = document_term_matrix, ...) %>%
           LDATS::LDA_select(...) %>%
           LDATS::MTS_prep(document_covariate_matrix) %>%
           LDATS::MTS_set(formula, nchangepoints, wts, ...) 
  return(out)
}

#' @title Calculate document weights (max value = 1)
#'
#' @param document_term_matrix matrix of documents (rows) by terms (columns)
#' @return vector of weights, one for each document, with the largest sample
#'   receiving a weight of 1.0
#'
#' @export
#'
doc_weights <- function(document_term_matrix){
  sample_sizes <- apply(document_term_matrix, 1, sum)
  out <- round(sample_sizes/max(sample_sizes), 3)  
  return(out)
}






