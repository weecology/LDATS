#' @title Run a set of Latent Dirichlet Allocation models
#' 
#' @description For a given dataset (counts of words across several 
#'   documents), conduct multiple Latent Dirichlet Allocation (LDA) models
#'   (using the Variational Expectation Maximization (VEM) algorithm) to 
#'   account for [1] uncertainty in the number of latent topics and
#'   [2] the impact of intial values in the estimation procedure 
#'   
#' @param MV Matrix of observation count data (rows: documents (\code{M}), 
#'   columns: terms (\code{V})).
#'
#' @param topics Vector of the number of topics to evaluate.
#'
#' @param nseeds Integer number of seeds (replicate starts) to use for each 
#'   value of \code{ntopics}.
#'
#' @param control Named list of control parameters to be used in 
#'   \code{\link[topicmodels]{LDA}} (note that "seed" will be overwritten).
#' 
#' @return List (class: "\code{LDA_list}") of LDA models (class: 
#'   "\code{LDA}").
#' 
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
#'   r_LDA <- LDA_set(MV = lda_data, topics = 2, nseeds = 2)                         
#' }
#' 
#' @export
#'
LDA_set <- function(MV, topics = 2, nseeds = 1, control = NULL){

  check_MV(MV)
  check_topics(topics)
  mod_topics <- rep(topics, each = length(seq(2, nseeds * 2, 2)))
  mod_seeds <- rep(seq(2, nseeds * 2, 2), length(topics))

  nmods <- length(mod_topics)
  mods <- vector("list", length = nmods)
  for (i in 1:nmods){
    control <- prep_LDA_control(seed = mod_seeds[i], control = control)
    mods[[i]] <- LDA(x = MV, k = mod_topics[i], control = control)
  }
  names(mods) <- paste0("c: ", mod_topics, ", seed: ", mod_seeds)
  class(mods) <- c("LDA_list", "list")  
  return(mods)
}

#' @title Verify that MV matrix is proper
#' 
#' @description Verify that the matrix of observations is conformable to
#'   a matrix of integers
#'   
#' @param MV Matrix of observation count data (rows: documents (\code{M}), 
#'   columns: terms (\code{V})).
#'
#' @return Nothing.
#' 
#' @export
#'
check_MV <- function(MV){
  MVm <- as.matrix(MV)
  if(!is.integer(MVm[1, 1])){
    stop("MV must be conformable to a matrix of integers")
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
  if (any(topics %% 1 != 0)){
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
#' @return List (class: "\code{LDAcontrol}") of LDA controls to be used.
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





