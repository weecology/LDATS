#' @title Determine the mode of a distribution
#'
#' @description Find the most common entry in a vector. Ties are not allowed,
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

#' @title Calculate document weights for a corpus
#'
#' @description Simple calculation of document weights based on the average
#'   number of words in a document within the corpus (mean value = 1).
#'
#' @param document_term_table Table of observation count data (rows: 
#'   documents, columns: terms. May be a class \code{matrix} or 
#'   \code{data.frame} but must be conformable to a matrix of integers,
#'   as verified by \code{\link{check_document_term_table}}.   
#'
#' @return Vector of weights, one for each document, with the average sample
#'   receiving a weight of 1.0.
#'
#' @export
#'
document_weights <- function(document_term_table){
  check_document_term_table(document_term_table)
  sample_sizes <- apply(document_term_table, 1, sum)
  round(sample_sizes/mean(sample_sizes), 3)  
}

#' @title Print with quieting
#'
#' @description Print a message (via \code{\link{cat}}) wrapped as in 
#'   \code{<wrapper><msg><wrapper>}, if desired.
#'
#' @param msg Message to be printed. \code{character}-class element.
#'
#' @param wrapper Wrapper \code{character} to use.
#'
#' @param quiet \code{logical} indicator of whether the message should be 
#'   printed.
#'
#' @export
#'
qprint <- function(msg, wrapper, quiet){
  if (quiet){
    return()
  }
  cat(paste0(wrapper, msg, wrapper, "\n"))
}


#' @title Create a properly symmetric variance covariance matrix
#'
#' @description A wrapper on \code{\link[stats]{vcov}} to produce a symmetric
#'   matrix. If the defauly matrix returned by \code{\link[stats]{vcov}} is
#'   symmetric it is returned simply. If it is not, in fact, symmetric
#'   (as occurs occasionally with \code{\link[nnet]{multinom}} applied to 
#'   proportions), the matrix is made symmetric by mirroring the lower
#'   triangle.
#'
#' @param x Model object that has a defined method for 
#'   \code{\link[stats]{vcov}}.
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
#' @param memoise_tf \code{logical} value indicatiing if \code{fun} should be 
#'   memoised.
#'
#' @return \code{fun}, memoised if desired.
#'
#' @export
#'
memoise_fun <- function(fun, memoise_tf){
  if (!("function" %in% class(fun))){
    stop("fun is not a function")
  }
  if (!("logical" %in% class(memoise_tf))){
    stop("memoise_tf is not logical")
  }
  if (memoise_tf){
    fun <- memoise(fun)
  }
  fun
}


#' @title Check that a control list is proper
#' 
#' @description Check that a list of controls is of the right class.
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


#' @title Check that document term table is proper
#' 
#' @description Check that the table of observations is conformable to
#'   a matrix of integers.
#'   
#' @param document_term_table Table of observation count data (rows: 
#'   documents, columns: terms. May be a class \code{matrix} or 
#'   \code{data.frame} but must be conformable to a matrix of integers,
#'   as verified by \code{\link{check_document_term_table}}. 
#' 
#' @export
#'
check_document_term_table <- function(document_term_table){
  document_term_table_m <- as.matrix(document_term_table)
  if(!is.integer(document_term_table_m[1, 1])){
    dtt <- "document_term_table"
    msg <- paste0(dtt, "is not conformable to a matrix of integers")
    stop(msg)
  }
}

#' @title Check that topics vector is proper
#' 
#' @description Check that the vector of numbers of topics is conformable to
#'   integers greater than 1.
#'   
#' @param topics Vector of the number of topics to evaluate for each model.
#'   Must be conformable to \code{integer} values.
#'
#' @export
#'
check_topics <- function(topics){
  if (!is.numeric(topics) || any(topics %% 1 != 0)){
    stop("topics vector must be integers")
  }
  if (any(topics < 2)){
    stop("minimum number of topics currently allowed is 2")
  }
}

#' @title Check that nseeds value or seeds vector is proper
#' 
#' @description Check that the vector of numbers of seeds is conformable to
#'   integers greater than 0.
#'   
#' @param nseeds \code{integer} number of seeds (replicate starts) to use for 
#'   each value of \code{topics} in the LDAs. Must be conformable to a  
#'   positive \code{integer} value.
#' 
#' @export
#'
check_seeds <- function(seeds){
  if (!is.numeric(seeds) || any(seeds %% 1 != 0)){
    stop("topics vector must be integers")
  }
}

