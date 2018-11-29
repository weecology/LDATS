#' @title Determine the mode of a distribution
#'
#' @description Find the most common entry in a vector. Ties are not included,
#'   the first value encountered within the modal set if there are ties is
#'   deemed the mode.
#'
#' @param x \code{numeric} vector.
#'
#' @return Numeric value of the mode.
#' 
#' @export 
#'
modalvalue <- function(x){
  as.numeric(names(sort(table(x), decreasing = TRUE))[1])
}

#' @title Calculate document weights on maximum number of words
#'
#' @description Simple calculation of document weights based on the maximum 
#'   number of words in a document within the corpus (max value = 1).
#'
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})). May be a class 
#'   \code{matrix} or \code{data.frame} but must be conformable to
#'   a code of integers. This table is a document-level summary of the data 
#'   noted as \eqn{w} (the word-level topic identity) in the math description. 
#'
#' @return Vector of weights, one for each document, with the largest sample
#'   receiving a weight of 1.0.
#'
#' @export
#'
document_weights <- function(document_term_table){
  sample_sizes <- apply(document_term_table, 1, sum)
  round(sample_sizes/max(sample_sizes), 3)  
}

#' @title Print with quieting
#'
#' @description Print a message (via \code{cat}) wrapped as in 
#'   \code{<wrapper><msg><wrapper>}.
#'
#' @param msg Message to be printed. \code{character}-class element.
#'
#' @param wrapper Wrapper \code{character} to use.
#'
#' @param quiet \code{logical} indicator of whether the message should be 
#'   printed.
#'
#' @return Nothing (message is printed, not returned).
#'
#' @export
#'
qprint <- function(msg, wrapper, quiet){
  if (quiet){
    return()
  }
  cat(paste0(wrapper, msg, wrapper, "\n"))
}