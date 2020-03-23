#' @title Check that arguments are properly formatted for usage
#'
#' @description Verify the class, structure, and values of inputted arguments
#'   to ensure proper LDATS modeling. \cr \cr
#'   \code{check_class} is a general class-verifier. \cr \cr
#'   \code{check_nonneg_integer} is a specified checking function for
#'     \code{numeric} values that must be integer-conformable and non-negative
#'     (0 is acceptable). \cr \cr
#'   \code{check_nonneg_integer_matrix} is a specified checking function for
#'     tables of \code{numeric} values that must be integer-conformable and 
#'     non-negative (0 is acceptable) and that must be conformable to 
#'     a \code{matrix}. \cr \cr
#'   \code{check_pos_integer} is a specified checking function for
#'     \code{numeric} values that must be integer-conformable and positive
#'     (0 is not acceptable). \cr \cr
#'   \code{check_control} verifies that a control \code{list} is a \code{list}
#'     \cr \cr
#'   \code{check_topics} ensures that the vector of numbers of topics is 
#'     positive integer-conformable. \cr \cr
#'   \code{check_replicates} ensures that the number of replicates is 
#'     positive integer-conformable. \cr \cr
#'   \code{check_nchangepoints} ensures that the of change points is 
#'     positive integer-conformable. \cr \cr
#'   \code{check_document_term_table} ensures that the table of document
#'     term counts is conformable to a \code{matrix} of positive integers. 
#'     \cr \cr
#'   \code{check_LDAs} verifies that the argument is either a
#'     \code{LDA} or \code{LDA_set} \code{list}. \cr \cr
#'   \code{check_document_covariate_table} check that the table of
#'    document-level covariates in the \code{LDAs} data is 
#'    conformable to a data frame and of the right size (correct number of 
#'    documents) for the document-topic output from the LDA models. \cr \cr
#'   \code{check_weights} ensures that the vector of document weights is 
#'    \code{numeric} and positive and inform the user if the average weight
#'    isn't 1.\cr \cr
#'   \code{check_timename} checks that the vector of time values is included 
#'     in the \code{document_covariate_table} and that it is either a 
#'     \code{integer}-conformable or a \code{Date}. 
#'     If it is a \code{Date}, the input is converted to an 
#'     \code{integer}, resulting in the timestep being 1 day, which is often
#'     not desired behavior. \cr \cr
#'   \code{check_formulas} verifies that the input contains only 
#'   \code{\link[stats]{formula}}s and that the response and predictor 
#'     variables are all included in \code{LDAs} data sets. \cr \cr
#'   \code{check_nchangepoints} checks that the \code{vector} of numbers of 
#'     changepoints is conformable to non-negative \code{integers}.\cr \cr
#'
#' @details
#' 
#' @param object An object whose class should be checked against 
#'   \code{eclass}.
#'
#' @param eclass Expected class of \code{object} to be checked. If more
#'   than one option is included, any are sufficient (\code{object} only
#'   needs to be one \code{eclass}, not all). 
#'   
#' @param control Control \code{list} to evaluate.
#'
#' @param topics \code{vector} of the number of topics to evaluate for each 
#'   model. Must be conformable to positive \code{integer} values.
#'
#' @param replicates \code{integer} number of replicate starts to use for 
#'   each value of \code{topics} in the LDAs. Must be conformable to  
#'   positive \code{integer} values.
#'
#' @param document_term_table Table of observation count data (rows: 
#'   documents, columns: terms. May be a \code{matrix} or 
#'   \code{data.frame} but must be conformable to a matrix of non-negative
#'   \code{integers}.
#'
#' @param nchangepoints \code{integer}-conformable \code{vector} of the 
#'   number of changepoints to evaluate (must be non-negative).
#'
#' @param weights \code{numeric} \code{vector} of the document weights to 
#'   evaluate, or \code{TRUE} for triggering internal weighting by document
#'   sizes.
#'
#' @param LDAs \code{LDA_models} \code{list} of LDA models or singular LDA 
#'   model (\code{LDA}) to evaluate.
#'
#' @param timename Column name for the time variable to evaluate in the
#'   \code{document_covariate_table} if provided.
#'
#' @param formulas \code{vector} of the \code{\link[stats]{formula}}s 
#'   to evaluate.
#'
#' @return an error message is thrown if the input is improper, otherwise 
#'   \code{NULL}.
#'
#'
#' @name argument_checking
#'



#' @rdname argument_checking
#'
#' @export
#'
check_class <- function(object, eclass = "list"){
  if(any(eclass == "nonneg_integer")){
    check_nonneg_integer(object)
    return(invisible())
  }
  if(any(eclass == "nonneg_integer_matrix")){
    check_nonneg_integer_matrix(object)
    return(invisible())
  }
  if(any(eclass == "pos_integer")){
    check_pos_integer(object)
    return(invisible())
  }
  if (!any(eclass %in% class(object))){
    object_name <- deparse(substitute(object))
    failed_eclass <- eclass[which(!(eclass %in% class(object)))]
    failed_eclass <- paste0(failed_eclass, collapse = " or ")
    stop(paste0(object_name, " is not a ", failed_eclass))
  }
}

#' @rdname argument_checking
#'
#' @export
#'
check_nonneg_integer <- function(object){
  object_name <- deparse(substitute(object))

  if (!is.numeric(object) || any(object %% 1 != 0)){
    stop(paste0(object_name, " must be integer-valued"))
  }
  if (any(object < 0)){
    stop(paste0(object_name, " must be non-negative"))
  }
}


#' @rdname argument_checking
#'
#' @export
#'
check_nonneg_integer_matrix <- function(object){
  object_name <- deparse(substitute(object))
  object_m <- as.matrix(object)
  check_nonneg_integer(object_m)
}



#' @rdname argument_checking
#'
#' @export
#'
check_pos_integer <- function(object){
  object_name <- deparse(substitute(object))
  if (!is.numeric(object) || any(object %% 1 != 0)){
    stop(paste0(object_name, " must be integer-valued"))
  }
  if (any(object <= 0)){
    stop(paste0(object_name, " must be positive"))
  }
}

#' @rdname argument_checking
#'
#' @export
#'
check_control <- function(control, eclass = "list"){
  check_class(object = control, eclass = eclass)
}


#' @rdname argument_checking
#'
#' @export
#'
check_topics <- function(topics){
  check_class(object = topics, eclass = "pos_integer")
}


#' @rdname argument_checking
#'
#' @export
#'
check_replicates <- function(replicates){
  check_class(object = replicates, eclass = "pos_integer")
}

#' @rdname argument_checking
#'
#' @export
#'
check_nchangepoints <- function(nchangepoints){
  check_class(object = nchangepoints, eclass = "nonneg_integer")
}

#' @rdname argument_checking
#'
#' @export
#'
check_document_term_table <- function(document_term_table){
  check_class(object = document_term_table, eclass = "nonneg_integer_matrix")

}


#' @rdname argument_checking
#'
#' @export
#'
check_weights <- function(weights){
  if(is.logical(weights)){
    if(weights){
      return(invisible())
    } else{
      stop("if logical, weights need to be TRUE")
    }   
  }
  if(!is.null(weights)){
    if (!is.numeric(weights)){
      stop("weights vector must be numeric")
    }
    if (any(weights <= 0)){
      stop("weights must be positive")
    }
    if (round(mean(weights)) != 1){
      warning("weights should have a mean of 1, fit may be unstable")
    }
  }
}

#' @rdname argument_checking
#'
#' @export
#'
check_LDAs <- function(LDAs){
  check_class(object = LDAs, eclass = c("LDA_set", "LDA"))
}




#' @rdname argument_checking
#'
#' @export
#'
check_timename <- function(LDAs, timename){
  
  document_cov_table <- LDAs[[1]][[1]]$data$train$document_covariate_table
  if (!("character" %in% class(timename))){
    stop("timename is not a character value")
  }
  if (length(timename) > 1){
    stop("timename can only be one value")
  }
  covariate_names <- colnames(document_cov_table)
  if ((timename %in% covariate_names) == FALSE){
    stop("timename not present in document covariate table")
  }
  time_covariate <- document_cov_table[ , timename]
  if (!(is.Date(time_covariate)) & 
      (!is.numeric(time_covariate) || !all(time_covariate %% 1 == 0))){
    stop("covariate indicated by timename is not an integer or a date")
  }
}






#' @rdname argument_checking
#'
#' @export
#'
check_formulas <- function(LDAs, formulas){

  dct <- LDAs[[1]][[1]]$data$train$document_covariate_table
  control <- LDAs[[1]][[1]]$control

  # response <- control$response

  if (!is(formulas, "list")) {
    if (is(formulas, "formula")) {
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain formula(s)")
    }
  } else if (!all(vapply(formulas, is, TRUE, "formula"))) {
      stop("formulas does not contain all formula(s)")
  }
  resp <- unlist(lapply(lapply(formulas, terms), attr, "response"))
  pred <- unlist(lapply(lapply(formulas, terms), attr, "term.labels"))
  if (any(resp != 0)) {
    stop("formula inputs should not include response variable")
  }
  if (!all(pred %in% colnames(dct))) {
    misses <- pred[which(pred %in% colnames(dct) == FALSE)]
    mis <- paste(misses, collapse = ", ")
    stop(paste0("formulas include predictors not present in data: ", mis))
  }
}





#' @rdname argument_checking
#'
#' @export
#'
check_document_covariate_table <- function(LDAs){

  dct <- LDAs[[1]][[1]]$data$train$document_covariate_table
  dtt <- LDAs[[1]][[1]]$data$train$document_term_table



  dct_df <- tryCatch(data.frame(dct),
                     warning = function(x){NA}, error = function(x){NA})
  if(is(LDAs, "LDA")){
    LDAs <- c(LDAs)
    class(LDAs) <- c("LDA_set", "list")
  }
  if (length(dct_df) == 1 && is.na(dct_df)){
    stop("document_covariate_table is not conformable to a data frame")
  }
  if (!is.null(LDAs)){
    if (nrow(data.frame(dct)) != nrow(LDAs[[1]][[1]]$document_topic_matrix)){
      stop("number of documents in covariate table is not equal to number of 
        documents observed")
    }
  } else if (!is.null(dtt)){
    if (nrow(data.frame(dct)) != nrow(data.frame(dtt))){
      stop("number of documents in covariate table is not equal to number of 
        documents observed")
    }
  }
}


