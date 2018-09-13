#' @title Run a set of Latent Dirichlet Allocation models
#' 
#' @description For a given dataset (counts of words across several 
#'   documents), conduct multiple Latent Dirichlet Allocation (LDA) models
#'   (using the Variational Expectation Maximization (VEM) algorithm; Blei et 
#'   al.) to account for [1] uncertainty in the number of latent topics and
#'   [2] the impact of intial values in the estimation procedure. 
#'
#' This function is a wrapper of the \code{\link[topicmodels]{LDA}} function
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
#' @param control Named list of control parameters to be used in 
#'   \code{\link[topicmodels]{LDA}} (note that "seed" will be overwritten).
#' 
#' @return List (class: "\code{LDA_list}") of LDA models (class: 
#'   "\code{LDA}").
#' 
#' @references 
#'   Blei, D. M., A. Y. Ng, and M. I. Jordan. 2003. Latent Dirichlet
#'   Allocation. Journal of Machine Learning Research 3:993-1022.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. Journal of Statistical Software 40:13.
#'
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
#'   r_LDA <- LDA_set(lda_data, topics = 2, nseeds = 2)                         
#' }
#' 
#' @export
#'
LDA_set <- function(document_term_table, topics = 2, nseeds = 1, 
                    control = NULL){

  check_document_term_table(document_term_table)
  check_topics(topics)
  mod_topics <- rep(topics, each = length(seq(2, nseeds * 2, 2)))
  mod_seeds <- rep(seq(2, nseeds * 2, 2), length(topics))
  dtt <- document_term_table

  nmods <- length(mod_topics)
  mods <- vector("list", length = nmods)
  for (i in 1:nmods){
    control <- prep_LDA_control(seed = mod_seeds[i], control = control)
    mods[[i]] <- LDA(x = dtt, k = mod_topics[i], control = control)
  }
  names(mods) <- paste0("c: ", mod_topics, ", seed: ", mod_seeds)
  class(mods) <- c("LDA_list", "list")  
  return(mods)
}

#' @title Verify that document term table is proper
#' 
#' @description Verify that the table of observations is conformable to
#'   a matrix of integers
#'   
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})).
#'
#' @return Nothing.
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
#' @param topics Vector of the number of topics to evaluate.
#'
#' @return Nothing.
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

#' @title Set the control inputs to include the seed
#' 
#' @description Update the control liste for the LDA model with the specific
#'   seed as indicated.
#'   
#' @param seed Integer number of seeds (replicate starts) to use for the 
#'   specific model.
#'
#' @param control Named list of control parameters to be used in 
#'   \code{\link[topicmodels]{LDA}} (note that "seed" will be overwritten).
#'
#' @return List (class: "\code{list}") of controls to be used in the LDA. 
#' 
#' @export
#'
prep_LDA_control <- function(seed, control = NULL){
  if(!is.list(control) & !is.null(control)){
    stop("control must be either a list or NULL")
  }
  if(is.list(control)){
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
#'   \code{LDA_list} object, based on a set of user-provided functions. The
#'   functions default to choosing the model with the lowest AIC value.
#'
#' @param LDA_models An object of class \code{LDA_list} produced by
#'   \code{LDA_set}.
#'
#' @param measurer,selector Function names for use in evaluation of the LDA
#'   models. \code{measurer} is used to create a value for each model
#'   and \code{selector} operates on the values to choose the model(s) to 
#'   pass on. 
#'
#' @return A reduced version of \code{LDA_models} that only includes the 
#'   selected LDA model(s). The returned object is still an object of
#'   class \code{LDA_list}.
#'
#' @examples
#' \dontrun{
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
#'   r_LDA <- LDA_set(lda_data, topics = 2, nseeds = 2)  
#'   select_LDA(r_LDA)                       
#' }
#'
#' @export
#'
select_LDA <- function(LDA_models = NULL, measurer = AIC, selector = min){

  if("LDA_list" %in% attr(LDA_models, "class") == FALSE){
    stop("LDA_models must be of class LDA_list")
  }
  
  lda_measured <- sapply(LDA_models, measurer) %>%
                  matrix(ncol = 1)
  lda_selected <- apply(lda_measured, 2, selector) 
  which_selected <- which(lda_measured %in% lda_selected)
  out <- lda_models[which_selected]
  class(out)  <- c("LDA_list", "list") 
  out
}
