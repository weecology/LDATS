soft_call <- function(fun = function(x){invisible(NULL)}, args = list(NULL), 
                      soften = FALSE){
  if(list_depth(args) == 0){
    args <- list(args)
  }
  if(soften){
    tryCatch(do.call(what = fun, args = args), 
             warning = function(x){eval(x$call)}, 
             error = function(x = list()){list(error = x$message)})
  } else{
    do.call(what = fun, args = args)
  }
}




time_order_data <- function(x, timename = "time"){
  time_order <- order(x[ , timename])
  x[time_order , ]
}


#' @title Initialize and tick through the progress bar
#'
#' @description \code{prep_pbar} creates and \code{update_pbar} steps
#'   through the progress bars (if desired) in, e.g., \code{\link{TS}}
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   iterative model.
#'
#' @param bar_type \code{character} value of possible types of progress bars.
#'   Currently available options are "rho" (for change point locations) and
#'   "eta" (for time series regressors).
#'
#' @param nr \code{integer} number of unique realizations, needed when
#'   \code{bar_type} = "eta".
#'
#' @param pbar The progress bar object returned from \code{prep_pbar}.
#'
#' @return 
#'   \code{prep_pbar}: the initialized progress bar object. \cr \cr
#'   \code{update_pbar}: the ticked-forward \code{pbar}.
#'
#' @export
#'
prep_pbar <- function(control = list(), bar_type = "rho", nr = NULL){
  if (!(bar_type %in% c("eta", "rho"))){
    stop("bar_type must be eta or rho")
  }
  if (!is.null(nr)){
    if (!is.numeric(nr) || any(nr %% 1 != 0)){
      stop("nr must be integer-valued")
    }
  }
  form <- "  [:bar] :percent eta: :eta"
  if (bar_type == "rho"){
    msg <- "    - estimating change point distribution"
    out <- progress_bar$new(form, control$nit, width = 60)
  }
  if (bar_type == "eta"){
    msg <- "    - estimating regressor distribution"
    out <- progress_bar$new(form, nr, width = 60)
  }
  messageq(msg, control$quiet)
  out
}

#' @rdname prep_pbar
#'
#' @export
#'
update_pbar <- function(pbar, control = list()){
  if (!("progress_bar" %in% class(pbar))){
    stop("pbar must be of class progress_bar")
  }
  if (control$quiet){
    return()
  }
  pbar$tick()
}


#' @title Determine the depth of a list
#'
#' @description Evaluate an input for the depth of its nesting. 
#'
#' @details If \code{xlist = list()}, then technically the input value is a 
#'  list, but is empty (of length \code{0}), so depth is returned as \code{0}.
#'
#' @param xlist Focal input \code{list}.
#'
#' @return \code{integer} value of the depth of the list.
#' 
#' @examples
#'  list_depth("a")
#'  list_depth(list())
#'  list_depth(list("a"))
#'  list_depth(list(list("a")))
#'
#' @export 
#'
list_depth <- function(xlist){
  xx <- match.call()
  xxx <- deparse(xx[[2]])
  if(xxx == "list()"){
    0L
  } else if (inherits(xlist, "data.frame")){
    0L
  } else if (is.list(xlist)){
    1L + max(sapply(xlist, list_depth))
  } else {
    0L
  }
}


#' @title Update a list's elements
#'
#' @description Update a list with new values for elements
#'
#' @param orig_list \code{list} to be updated with \code{...}. 
#'
#' @param ... Named elements to update in \code{orig_list}
#'
#' @return Updated \code{list}.
#'
#' @examples
#'  orig_list <- list(a = 1, b = 3, c = 4)
#'  update_list(orig_list)
#'  update_list(orig_list, a = "a")
#'  update_list(orig_list, a = 10, b = NULL)
#'
#' @export
#'
update_list <- function(orig_list = list(), ...){
  if(!is.list(orig_list)){
    stop("orig_list must be a list", call. = FALSE)
  } 
  update_elems <- list(...)
  nupdate_elems <- length(update_elems)
  norig_elems <- length(orig_list)
  update_list <- vector("list", length = norig_elems)
  names(update_list) <- names(orig_list)
  if(norig_elems > 0){
    for(i in 1:norig_elems){
      if(!is.null(orig_list[[i]])){
        update_list[[i]] <- orig_list[[i]]
      }
    }
  }
  if(nupdate_elems > 0){
    names_update_elems <- names(update_elems)
    for(i in 1:nupdate_elems){
      if(!is.null(update_elems[[i]])){
        update_list[[names_update_elems[i]]] <- update_elems[[i]]
      }
    }
  }
  update_list
}



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
  if (!(is.null(memoise_tf) || "logical" %in% class(memoise_tf))){
    stop("memoise_tf is not logical")
  }
  if (!is.null(memoise_tf) && memoise_tf){
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

#' @title Check that nseeds value or vector is proper
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
#'  check_nseeds(1)
#'  check_nseeds(2)
#'
#' @export
#'
check_nreps <- function(nreps){
  if (!is.numeric(nreps) || any(nreps %% 1 != 0)){
    stop("nreps must be integer-conformable")
  }
  return()
}

# provides a functionality that can be used in testing for non-symmetric
# vcov matrix
vcov.dummy <- function(object, ...){
  matrix(c(1, 2, 2.1, 3), 2, 2)
}
