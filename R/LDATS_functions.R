
#' @title Two-stage LDA-time series analysis
#'
#' @param data data set (currently just the data for the LDA)
#' @param LDA_eval function name for evaluation of the LDA models
#' @param LDA_selector function name for selecting the LDA models 
#' @param ... additional arguments to be passed to subfunctions
#' @return (currently) the output of the selected LDA model
#'
#' @export
#'
LDA_TS <- function(data = NULL, LDA_eval = quote(AIC),  
                   LDA_selector = quote(min), ...){
  lda_mods <- LDATS::LDA(data, ...)

  lda_eval <- sapply(lda_mods, LDA_eval) %>%
              matrix(ncol = 1)
  lda_selected <- apply(lda_eval, 2, LDA_selector) 
  which_selected <- which(lda_eval %in% lda_selected)
  selected_lda <- lda_mods[[which_selected]]
  out <- selected_lda
  return(out)
}