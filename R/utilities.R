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