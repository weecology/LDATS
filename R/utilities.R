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
  if (!is.numeric(x)){
    stop("x must be numeric")
  }
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
  check_document_term_table(document_term_table)
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


#' @title Create a properly symmetrix variance covariance matrix
#'
#' @description If the default matrix returned by \code{vcov} is not, in fact,
#'   symmetric (as occurs occasionally in the \code{multinom} function from
#'   \code{nnet}), make the matrix symmetric by mirroring the lower triangle.
#'
#' @param x Model object that has a defined method for \code{vcov}.
#'
#' @return Properly symmetric variance covariance matrix.
#'
#' @export
#'
mirror_vcov <- function(x){
  
  vcv <- tryCatch(vcov(x), error = function(x){NULL})
  if (is.null(vcv)){
    stop("`vcov` not defined for x")
  }
  if (isSymmetric(vcv)){
    return(vcv)
  }
  lt <- vcv[lower.tri(vcv)]
  ut <- t(vcv)[lower.tri(t(vcv))]
  if (any((abs(lt - ut) > (0.001 * abs(lt))))){
    stop("Relative discrepancies in model vcov matrix are too large.")  
  }
  vcv2 <- t(vcv)
  vcv2[lower.tri(vcv2)] <- vcv[lower.tri(vcv)]
  vcv2
}

#' @title Normalize a vector
#'
#' @description Normalize a \code{numeric} vector to be on the scale of [0,1].
#'
#' @param x \code{numeric} vector.
#'
#' @return Normalized \code{x}.
#' 
#' @export 
#'
normalize <- function(x){
  if (!is.numeric(x)){
    stop("x must be numeric")
  }
  (x - min(x)) / (max(x) - min(x))
}


#' @title Logical control on whether or not to memoise
#'
#' @description This function provides a simple, logical toggle control on
#'   whether the function \code{fun} should be memoised via
#'   \code{\link[memoise]{memoise}} or not. 
#'
#' @param fun Function name to (potentially) be memoised.
#'
#' @param memoise_fun \code{logical} value indicatiing if \code{fun} should be 
#'   memoised.
#'
#' @return \code{fun}, memoised if desired.
#'
#' @export
#'
memoise_fun <- function(fun, memoise_fun){
  if (!("function" %in% class(fun))){
    stop("fun is not a function")
  }
  if (!("logical" %in% class(memoise_fun))){
    stop("memoise_fun is not logical")
  }
  if (memoise_fun){
    fun <- memoise(fun)
  }
  fun
}


#' @title Verify that a control list is proper
#' 
#' @description Verify that a list of controls is of the right class.
#'   
#' @param control Control list to evaluate.
#' 
#' @param eclass Expected class of the list to be evaluated.
#'
#' @export
#'
check_control <- function(control, eclass = "TS_controls"){
  if (!(eclass %in% class(control))){
    stop(paste0("control is not of class ", eclass))
  }
}
