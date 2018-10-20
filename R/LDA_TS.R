#' @title Run a full set of Latent Dirichlet Allocation models then Time 
#'   Series models
#'
#' @description This function runs a complete suite of Latent Dirichlet
#'   Allocation (LDA) models (via \code{\link{LDA_set}}), selects the choice
#'   LDA model(s) (via \code{\link{select_LDA}}), runs a complete set of Time 
#'   Series (TS) models (via \code{\link{TS_on_LDA}}) on the choice LDA
#'   model(s), then selects the best TS model (via \code{select_TS}). 
#' 
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})). May be a class 
#'   \code{matrix} or \code{data.frame} but must be conformable to
#'   a code of integers. This table is a document-level summary of the data 
#'   noted as \eqn{w} (the word-level topic identity) in the math description. 
#'
#' @param topics Vector of the number of topics to evaluate.
#'
#' @param nseeds Integer number of seeds (replicate starts) to use for each 
#'   value of \code{topics}.
#'
#' @param lda_control Named list of control parameters to be used in 
#'   \code{\link[topicmodels]{LDA}} (note that "seed" will be overwritten).
#'
#' @return Nothing yet. 
#'
#' @export
#'
LDA_TS <- function(document_term_table, topics = 2, nseeds = 1, 
                   lda_control = NULL){


}