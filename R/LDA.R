#' @title Run a set of Latent Dirichlet Allocation models
#' 
#' @description For a given dataset (counts of words across several 
#'   documents), conduct multiple Latent Dirichlet Allocation (LDA) models
#'   (using the Variational Expectation Maximization (VEM) algorithm; Blei  
#'   \emph{et al.}) to account for [1] uncertainty in the number of latent 
#'   topics and [2] the impact of intial values in the estimation procedure. 
#'
#' This function is a wrapper of the \code{LDA} function
#'   in the \code{topicmodels} package (Grun and Hornik 2011).
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
#' @param control Class \code{LDA_controls} list of control parameters to be
#'   used in \code{LDA} (note that \code{seed} will be overwritten).
#' 
#' @return List (class: \code{LDA_set}) of LDA models (class: \code{LDA}).
#' 
#' @references 
#'   Blei, D. M., A. Y. Ng, and M. I. Jordan. 2003. Latent Dirichlet
#'   Allocation. \emph{Journal of Machine Learning Research} 
#'   \strong{3}:993-1022.
#'   \href{http://jmlr.csail.mit.edu/papers/v3/blei03a.html}{link}.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#'
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   r_LDA <- LDA_set(lda_data, topics = 2, nseeds = 2)                         
#' }
#' 
#' @export
#'
LDA_set <- function(document_term_table, topics = 2, nseeds = 1, 
                    control = LDA_controls_list()){
  check_LDA_set_inputs(document_term_table, topics, nseeds, control)
  mod_topics <- rep(topics, each = length(seq(2, nseeds * 2, 2)))
  mod_seeds <- rep(seq(2, nseeds * 2, 2), length(topics))
  nmods <- length(mod_topics)
  mods <- vector("list", length = nmods)
  for (i in 1:nmods){
    LDA_msg(mod_topics[i], mod_seeds[i], control)
    control_i <- prep_LDA_control(seed = mod_seeds[i], control = control)
    mods[[i]] <- LDA(document_term_table, k = mod_topics[i], 
                     control = control_i)
  }
  package_LDA_set(mods, mod_topics, mod_seeds)
}

#' @title Calculate the log likelihood of a VEM LDA model fit
#'
#' @description Imported calculations from topicmodels, as applied to
#'   Latent Dirichlet Allocation fit with Variational Expectation Maximization
#'   via \code{LDA}. 
#'
#' @param object A \code{LDA_VEM}-class object.
#'
#' @param ... Not used, simply included to maintain method compatability.
#'
#' @return Log likelihood of the model, also with \code{df} (degrees of
#'    freedom) and \code{nobs} (number of observations) values.
#'
#' @references 
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#'
#' @export
#'
logLik.LDA_VEM <- function(object, ...){
  val <- sum(object@loglikelihood)
  df <- as.integer(object@control@estimate.alpha) + length(object@beta)
  attr(val, "df") <- df
  attr(val, "nobs") <- object@Dim[1]
  class(val) <- "logLik"
  val
}

#' @title Verify that all of the inputs are proper for LDA_set 
#' 
#' @description Verify that the table of observations is conformable to
#'   a matrix of integers, the number of topics is an integer, the number of 
#'   seeds is an integer and the controls list is proper.
#'   
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})).
#'
#' @param topics Vector of the number of topics (\eqn{k}) to evaluate.
#'
#' @param nseeds Integer number of seeds (replicate starts) to use for each 
#'   value of \code{topics}.
#'
#' @param control Class \code{LDA_controls} list of control parameters to be
#'   used in \code{LDA} (note that \code{seed} will be overwritten).
#' 
#' @export
#'
check_LDA_set_inputs <- function(document_term_table, topics, nseeds, 
                                 control){
  check_document_term_table(document_term_table)
  check_topics(topics)
  check_seeds(nseeds)
  check_control(control, "LDA_controls")
}

#' @title Verify that document term table is proper
#' 
#' @description Verify that the table of observations is conformable to
#'   a matrix of integers
#'   
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})).
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

#' @title Verify that topics vector is proper
#' 
#' @description Verify that the vector of numbers of topics is conformable to
#'   integers greater than 1.
#'   
#' @param topics Vector of the number of topics to evaluate (\eqn{k}).
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

#' @title Verify that nseeds value or seeds vector is proper
#' 
#' @description Verify that the vector of numbers of seeds is conformable to
#'   integers greater than 1.
#'   
#' @param seeds Value of the number of random seeds to evaluate.
#'
#' @return Nothing.
#' 
#' @export
#'
check_seeds <- function(seeds){
  if (!is.numeric(seeds) || any(seeds %% 1 != 0)){
    stop("topics vector must be integers")
  }
}

#' @title Set the control inputs to include the seed
#' 
#' @description Update the control liste for the LDA model with the specific
#'   seed as indicated.
#'   
#' @param seed Integer number of seeds (replicate starts) to use for the 
#'   specific model.
#'
#' @param control Named list of control parameters to be used in 
#'   \code{LDA} (note that "seed" will be overwritten).
#'
#' @return List (class: "\code{list}") of controls to be used in the LDA. 
#' 
#' @export
#'
prep_LDA_control <- function(seed, control = NULL){
  if("LDA_controls" %in% class(control)){
    class(control) <- "list"
    control$quiet <- NULL
    control$measurer <- NULL
    control$selector <- NULL
    control$seed <- seed
  }
  if(is.null(control)){
    control <- list(seed = seed)
  }
  control
}

#' @title Select the best LDA model(s) for use in time series
#'
#' @description Select the best model(s) of interest from an
#'   \code{LDA_set} object, based on a set of user-provided functions. The
#'   functions default to choosing the model with the lowest AIC value.
#'
#' @param LDA_models An object of class \code{LDA_set} produced by
#'   \code{LDA_set}.
#'
#' @param control Class \code{LDA_controls} list including (named) elements
#'   corresponding to the \code{measurer} and \code{evaluator} functions.
#'
#' @return A reduced version of \code{LDA_models} that only includes the 
#'   selected LDA model(s). The returned object is still an object of
#'   class \code{LDA_set}.
#'
#' @examples
#' \dontrun{
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   r_LDA <- LDA_set(lda_data, topics = 2, nseeds = 2)  
#'   select_LDA(r_LDA)                       
#' }
#'
#' @export
#'
select_LDA <- function(LDA_models = NULL, control = LDA_controls_list()){

  measurer <- control$measurer
  selector <- control$selector
  if("LDA_set" %in% attr(LDA_models, "class") == FALSE){
    stop("LDA_models must be of class LDA_set")
  }
  
  lda_measured <- sapply(LDA_models, measurer) %>%
                  matrix(ncol = 1)
  lda_selected <- apply(lda_measured, 2, selector) 
  which_selected <- which(lda_measured %in% lda_selected)
  out <- LDA_models[which_selected]
  class(out)  <- c("LDA_set", "list") 
  out
}

#' @title Package the output from LDA_set
#'
#' @description Name the elements (LDA models) and set the class 
#'   (\code{LDA_set}) of the models returned by the \code{LDA_set} function.
#'
#' @param mods Fitted models returned from \code{LDA}.
#'
#' @param mod_topics Vector of \code{integer} values corresponding to the 
#'   number of topics in each model.
#' 
#' @param mod_seeds Vector of \code{integer} values corresponding to the 
#'   seed used for each model.
#'
#' @return List (class: \code{LDA_set}) of LDA models (class: \code{LDA}).
#'
#' @export
#'
package_LDA_set <- function(mods, mod_topics, mod_seeds){
  if (!("LDA_VEM" %in% class(mods[[1]]))){
    stop("mods not of class LDA_VEM")
  }
  check_topics(mod_topics)
  if (!is.numeric(mod_seeds) || any(mod_seeds%% 1 != 0)){
    stop("mod_seeds must be integers")
  }
  names(mods) <- paste0("k: ", mod_topics, ", seed: ", mod_seeds)
  class(mods) <- c("LDA_set", "list")  
  mods
}

#' @title Create the model-running-message for an LDA
#'
#' @description Produce and print the message for a given LDA model.
#'
#' @param mod_topics \code{integer} value corresponding to the number of 
#'   topics in the model.
#' 
#' @param mod_seeds \code{integer} value corresponding to the seed used for 
#'   the model.
#'
#' @param control Class \code{LDA_controls} list of control parameters to be
#'   used in \code{LDA} (note that "seed" will be overwritten).
#'
#' @export
#'
LDA_msg <- function(mod_topics, mod_seeds, control){
  check_topics(mod_topics)
  check_seeds(mod_seeds)
  topic_msg <- paste0("Running LDA with ", mod_topics, " topics ")
  seed_msg <- paste0("(seed ", mod_seeds, ")")
  qprint(paste0(topic_msg, seed_msg), "", control$quiet)
}

#' @title Create control list for LDA model
#'
#  @description This function provides a simple creation and definition of the
#'   list used to control the set of LDA models. It is set up to be easy to 
#'   work with the existing control capacity of the \code{LDA} function.
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly.
#'
#' @param measurer,selector Function names for use in evaluation of the LDA
#'   models. \code{measurer} is used to create a value for each model
#'   and \code{selector} operates on the values to choose the model(s) to 
#'   pass on. 
#'
#' @param ... Additional arguments to be passed to \code{LDA} as a 
#'   \code{control} input.
#'
#' @return Class \code{LDA_controls} list for controlling the LDA model fit.
#'
#' @export
#'
LDA_controls_list <- function(quiet = FALSE, measurer = AIC, selector = min,
                              ...){
  out <- list(quiet = quiet, measurer = measurer, selector = selector, ...)
  class(out) <- c("LDA_controls", "list")
  out
}
