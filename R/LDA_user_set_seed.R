#' @title Run a set of Latent Dirichlet Allocation models with user-specified seed
#' 
#' @description For a given dataset consisting of counts of words across 
#'   multiple documents in a corpus, conduct multiple Latent Dirichlet 
#'   Allocation (LDA) models (using the Variational Expectation 
#'   Maximization (VEM) algorithm; Blei \emph{et al.} 2003) to account for [1]  
#'   uncertainty in the number of latent topics and [2] the impact of initial
#'   values in the estimation procedure. \cr \cr
#'   \code{LDA_set} is a list wrapper of \code{\link[topicmodels]{LDA}}
#'   in the \code{topicmodels} package (Grun and Hornik 2011). \cr \cr
#'   \code{check_LDA_set_inputs} checks that all of the inputs 
#'   are proper for \code{LDA_set} (that the table of observations is 
#'   conformable to a matrix of integers, the number of topics is an integer, 
#'   the number of seeds is an integer and the controls list is proper).
#'   
#' @param document_term_table Table of observation count data (rows: 
#'   documents, columns: terms. May be a class \code{matrix} or 
#'   \code{data.frame} but must be conformable to a matrix of integers,
#'   as verified by \code{\link{check_document_term_table}}.   
#'  
#' @param topics Vector of the number of topics to evaluate for each model.
#'   Must be conformable to \code{integer} values.
#'
#' @param seed Seed to use for each 
#'   value of \code{topics}. 
#'
#' @param control A \code{list} of parameters to control the running and 
#'   selecting of LDA models. Values not input assume default values set 
#'   by \code{\link{LDA_set_control}}. Values for running the LDAs replace 
#'   defaults in (\code{LDAcontol}, see \code{\link[topicmodels]{LDA}} (but if
#'    \code{seed} is given, it will be overwritten; use \code{iseed} instead).
#' 
#' @return 
#'   \code{LDA_set}: \code{list} (class: \code{LDA_set}) of LDA models 
#'   (class: \code{LDA_VEM}).
#'   \code{check_LDA_set_inputs}: an error message is thrown if any input is 
#'   improper, otherwise \code{NULL}.
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
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   r_LDA <- LDA_set_user_seeds(lda_data, topics = 2, seed = 2)                         
#' 
#' @export
#'
LDA_set_user_seeds <- function(document_term_table, topics = 2, seed = 1, 
                    control = list()){
  nseeds = length(seed)
  check_LDA_set_inputs(document_term_table, topics, nseeds = nseeds, control)
  control <- do.call("LDA_set_control", control)
  mod_topics <- rep(topics, each = length(seq(2, length(seed) * 2, 2)))
  iseed <- control$iseed
  #mod_seeds <- rep(seq(iseed, iseed + (nseeds - 1)* 2, 2), length(topics))
  mod_seeds <- seed
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