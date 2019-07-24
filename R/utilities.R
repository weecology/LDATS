#' @title Calculate the log-sum-exponential (LSE) of a vector
#'
#' @description Calculate the exponent of a vector (offset by the max), sum
#'   the elements, calculate the log, remove the offset. 
#'
#' @param x \code{numeric} vector
#' 
#' @return The LSE. 
#'
#' @examples
#'   logsumexp(1:10)
#'
#' @export
#'
logsumexp <- function(x){
  if(!is.numeric(x)){
    stop("x must be numeric")
  }
  y <- max(x)
  y + log(sum(exp(x - y)))
}

#' @title Calculate the softmax of a vector or matrix of values
#'
#' @description Calculate the softmax (normalized exponential) of a vector
#'  of values or a set of vectors stacked rowwise. 
#'
#' @param x \code{numeric} vector or matrix
#' 
#' @return The softmax of \code{x}. 
#'
#' @examples
#'   dat <- matrix(runif(100, -1, 1), 25, 4)
#'   softmax(dat)
#'   softmax(dat[,1])
#'
#' @export
#'
softmax <- function(x){
  if(!is.numeric(x)){
    stop("x must be numeric")
  }
  if(length(dim(x)) == 0){
    exp(x - logsumexp(x))
  } else if (length(dim(x)) == 2){
    outm <- matrix(NA, nrow = nrow(x), ncol = ncol(x))
    for(i in 1:nrow(x)){
      outm[i,] <- exp(x[i,] - logsumexp(x[i,]))
    }
    outm
  } else{
    stop("currently only defined for vector or matrix")
  }
}

#' @title Calculate AICc
#'
#' @description Calculate the small sample size correction of
#'   \code{\link{AIC}} for the input object. 
#'
#' @param object An object for which \code{\link{AIC}} and 
#'   \code{\link{logLik}} have defined methods.
#'
#' @return \code{numeric} value of AICc.
#' 
#' @examples
#'   dat <- data.frame(y = rnorm(50), x = rnorm(50))
#'   mod <- lm(dat)
#'   AICc(mod)
#'
#' @export 
#'
AICc <- function(object){
  aic <- AIC(object)
  ll <- logLik(object)
  np <- attr(ll, "df")
  no <- attr(ll, "nobs")
  aic + (2 * np^2 + 2 * np)/(no - np - 1)
}

#' @title Replace if TRUE
#'
#' @description If the focal input is \code{TRUE}, replace it with 
#'   alternative. 
#'
#' @param x Focal input.
#'
#' @param alt Alternative value.
#'
#' @return \code{x} if not \code{TRUE}, \code{alt} otherwise.
#' 
#' @examples
#'  iftrue()
#'  iftrue(TRUE, 1)
#'  iftrue(2, 1)
#'  iftrue(FALSE, 1)
#'
#' @export 
#'
iftrue <- function(x = TRUE, alt = NULL){
  if (is.logical(x) && x){
    x <- alt
  }
  x
}

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
#' @examples
#'  d1 <- c(1, 1, 1, 2, 2, 3)
#'  modalvalue(d1)
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
#' @examples
#'  data(rodents)
#'  document_weights(rodents$document_term_table)
#'
#' @export
#'
document_weights <- function(document_term_table){
  check_document_term_table(document_term_table)
  sample_sizes <- apply(document_term_table, 1, sum)
  round(sample_sizes/mean(sample_sizes), 3)  
}

#' @title Optionally generate a message based on a logical input
#'
#' @description Given the input to \code{quiet}, generate the message(s) 
#'   in \code{msg} or not.
#'
#' @param msg \code{character} vector of the message(s) to generate or 
#'   \code{NULL}. If more than one element is contained in \code{msg}, they
#'   are concatenated with a newline between.
#'
#' @param quiet \code{logical} indicator controlling if the message is
#'   generated.
#'
#' @examples
#'   messageq("hello")
#'   messageq("hello", TRUE)
#'
#' @export
#'
messageq <- function(msg = NULL, quiet = FALSE){
  if (!quiet){
    msg2 <- paste(msg, collapse = "\n")
    message(msg2)
  }
}


#' @title Create a properly symmetric variance covariance matrix
#'
#' @description A wrapper on \code{\link[stats]{vcov}} to produce a symmetric
#'   matrix. If the default matrix returned by \code{\link[stats]{vcov}} is
#'   symmetric it is returned simply. If it is not, in fact, symmetric
#'   (as occurs occasionally with \code{\link[nnet]{multinom}} applied to 
#'   proportions), the matrix is made symmetric by averaging the lower and
#'   upper triangles. If the relative difference between the upper and lower 
#'   triangles for any entry is more than 0.1% a warning is thrown. 
#'
#' @param x Model object that has a defined method for 
#'   \code{\link[stats]{vcov}}.
#'
#' @return Properly symmetric variance covariance \code{matrix}.
#'
#' @examples
#'   dat <- data.frame(y = rnorm(50), x = rnorm(50))
#'   mod <- lm(dat)
#'   mirror_vcov(mod)
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
  tvcv <- t(vcv)
  lt <- vcv[lower.tri(vcv)]
  ut <- tvcv[lower.tri(tvcv)]
  if (any((abs(lt - ut) > (0.001 * abs(lt))))){
    warning("Relative discrepancies in model vcov matrix are large.")  
  }
  bound <- cbind(as.numeric(vcv), as.numeric(tvcv))
  avg <- apply(bound, 1, mean)
  vcv2 <- matrix(avg, nrow(vcv), ncol(vcv))
  rownames(vcv2) <- rownames(vcv)
  colnames(vcv2) <- colnames(vcv)
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
#' @examples
#'  normalize(1:10)
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
#' @param memoise_tf \code{logical} value indicating if \code{fun} should be 
#'   memoised.
#'
#' @return \code{fun}, memoised if desired.
#'
#' @examples
#'   sum_memo <- memoise_fun(sum)
#'
#' @export
#'
memoise_fun <- function(fun, memoise_tf = TRUE){
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
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  check_control(list())
#'
#' @export
#'
check_control <- function(control, eclass = "list"){
  if (!(eclass %in% class(control))){
    stop(paste0("control is not a ", eclass))
  }
  return()
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
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  data(rodents)
#'  check_document_term_table(rodents$document_term_table)
#'
#' @export
#'
check_document_term_table <- function(document_term_table){
  document_term_table_m <- as.matrix(document_term_table)
  if(any(document_term_table_m %% 1 != 0)){
    dtt <- "document_term_table"
    msg <- paste0(dtt, " is not conformable to a matrix of integers")
    stop(msg)
  }
  return()
}

#' @title Check that topics vector is proper
#' 
#' @description Check that the vector of numbers of topics is conformable to
#'   integers greater than 1.
#'   
#' @param topics Vector of the number of topics to evaluate for each model.
#'   Must be conformable to \code{integer} values.
#'
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  check_topics(2)
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
  return()
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
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#' @examples
#'  check_seeds(1)
#'  check_seeds(2)
#'
#' @export
#'
check_seeds <- function(nseeds){
  if (!is.numeric(nseeds) || any(nseeds %% 1 != 0)){
    stop("nseeds vector must be integers")
  }
  return()
}

# provides a functionality that can be used in testing for non-symmetric
# vcov matrix
vcov.dummy <- function(object, ...){
  matrix(c(1, 2, 2.1, 3), 2, 2)
}
