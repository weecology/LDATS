#' @title Two-stage LDA-time series analysis
#'
#' @description Run LDA, select LDA, run MTS; return all outputs
#'
#' @param document_term_matrix matrix of documents (rows) by terms (columns)
#' @param document_covariate_matrix matrix of documents (rows) by covariates
#'   (columns)
#' @param formula formula or vector of formulas for the continuous change
#' @param nchangepoints vector of the number of change points to include in 
#'   the model
#' @param ... additional arguments to be passed to subfunctions
#' @return a list of [1] the LDA model(s), [2] the selected LDA model(s), and
#'   [3] the time series model(s) on the selected LDA model(s)
#'
#' @examples 
#'   \dontrun{
#'     data(rodents)
#'     lda_data <- select(rodents, -c(newmoon, date, plots, traps))
#'     ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])
#'     r_LDATS <- LDA_TS(lda_data, ts_data, formula = ~time, ntopics = 2:5,
#'                       nseeds = 2, ncores = 4, nit = 100) 
#'   }
#'
#' @export
#'
LDA_TS <- function(document_term_matrix = NULL, 
                   document_covariate_matrix = NULL, 
                   formula = ~1, nchangepoints = 1, ...){

  weights <- doc_weights(document_term_matrix)
  ldas <- LDA(data = document_term_matrix, ...) 
  selected <- LDA_select(ldas, ...) 
  mtss <- selected %>%
          MTS_prep(document_covariate_matrix) %>%
          MTS_set(formula, nchangepoints, weights, ...) 

  out <- list(ldas, selected, mtss)
  names(out) <- c("LDA model(s)", "Selected LDA model(s)", "MTS model(s)")
}

#' @title Calculate document weights
#'
#' @description Simple calculation based on the maximum (max value = 1)
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
