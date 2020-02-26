#' @title Run a set of Latent Dirichlet Allocation models
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
#'   Must be conformable to \code{integer} values. REDUCED DIMENSIONALITY
#'
#' @param nreps Number of replicate starts to use for each 
#'   value of \code{topics}. Must be conformable to \code{integer} value.
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
#'   r_LDA <- LDA(lda_data, topics = 2, nseeds = 2)                         
#' 
#' @export
#'
LDA <- function(data, topics = 2:3, reps_per_topic = 2, control = list()){
 # check_args()
  control <- do.call("LDA_control", control)

  messageq("----Linguistic Decomposition Analyses----", control$quiet)

  LDAs <- prep_LDA_models(data = data, topics = topics, 
                          reps_per_topic = reps_per_topic,
                          control = control)
  nLDA <- length(LDA_models)

  for (i in 1:nLDA){
    LDAs[[i]] <- LDA_call(data = LDAs[[i]]$data, topics = LDAs[[i]]$topics,
                          rep = LDAs[[i]]$rep, control = control)
  }
  # eval_LDA
  # select_LDA
  package_LDA(selected_LDAs = selected_LDAs, 
              LDAs = LDAs, LDA_topics = LDA_topics, LDA_seeds = LDA_seeds)
}

prep_LDA_models <- function(data, topics = 2:3, reps_per_topic = 2,
                            control = list()){
  data <- conform_data(data = data, control = control)
  subsets <- names(data)
  LDA_topics <- rep(topics, each = length(seq(2, reps_per_topic * 2, 2)))
  LDA_reps <- rep(seq(1, 1 + (reps_per_topic - 1), 1), length(topics))
  LDA_subsets <- rep(subsets, each = length(LDA_reps))  
  LDA_reps <- rep(LDA_reps, length(subsets))
  LDA_topics <- rep(LDA_topics, length(subsets))
  nLDA <- length(LDA_topics)
  LDAs <- vector("list", length = nLDA)
  for(i in 1:nLDA){
    LDAs[[i]] <- list(data = data[[LDA_subsets[[i]]]], 
                      topics = LDA_topics[[i]], rep = LDA_reps[[i]])
  }
  names(LDAs) <- paste0("model_", 1:nLDA)

  LDAs
}


select_LDA <- function(LDAs = list(), control = list()){
  vals <- measure_LDA(LDAs = LDAs, control = control)
  
}

#working in here
#coming along!
#but remember that we'll need to manage passing data around for test/train
# and the multiple versions of data

measure_LDA <- function(LDAs = list(), control = list()){
  fun <- control$measurer_function
  args <- control$measurer_args
  nLDAs <- length(LDAs)
  vals <- rep(NA, nLDAs)
  for(i in 1:nLDAs){
    args <- update_list(args, object = LDAs[[i]])
    vals[i] <- do.call(what = fun, args = args)
  }
  vals
}

LDA_call <- function(data, topics = 2, rep = 1, control = list()){
  control <- do.call("LDA_control", control)
  LDA_msg(topics = topics, rep = rep, quiet = control$quiet)
  fun <- control$LDA_function
  args <- update_list(control$LDA_args, data = data, topics = topics, 
                      rep = rep)
  if(control$soften){
    tryCatch(do.call(what = fun, args = args), 
             warning = LDA_call_soft_warning, 
             error = LDA_call_soft_error)
  } else{
    do.call(what = fun, args = args)
  }
}

AIC.LDA <- function(object, ..., k = 2){
  lls <- logLik(object)
  -2 * as.numeric(lls) + k * attr(lls, "df")
}

logLik.LDA <- function(object, ...){
  object$log_likelihood
}

LDA_call_soft_warning <- function(x){
  eval(x$call)
}

LDA_call_soft_error <- function(x = list()){
  list(params = list(), document_topic_matrix = data.frame(), 
       log_likelihood = NA, error = x$message)
}


topicmodels_LDA <- function(data, topics = 2, rep = 1, method = "VEM", 
                            seeded = TRUE, ...){
  if(topics == 1){
    identity_LDA(data)
  } else{
    fun_control <- list(...)
    if(seeded){
      fun_control <- update_list(fun_control, seed = rep * 2)
    }
    mod <- topicmodels::LDA(x = data$document_term_table, k = topics, 
                            method = method, control = fun_control)
    mod_ll <- sum(mod@loglikelihood)
    df <- as.integer(mod@control@estimate.alpha) + length(mod@beta)
    attr(mod_ll, "df") <- df
    attr(mod_ll, "nobs") <- mod@Dim[1] * mod@Dim[2]
    class(mod_ll) <- "logLik"
    out <- list(params = list(alpha = mod@alpha, beta = mod@beta),
                document_topic_matrix = mod@gamma, 
                log_likelihood = mod_ll, data = data)
    class(out) <- c("LDA", "list")
    out
  } 
}

identity_LDA <- function(data, topics = NULL, rep = 1){
  nterms <- NCOL(data$document_term_table)
  if(is.null(topics)){
    nterms <- topics
  } else if(topics != nterms){
    stop("number of topics does not equal number of terms")
  }
  document_topic_table <- data$document_term_table 
  document_topic_table <- document_topic_table / rowSums(document_topic_table)
  colnames(document_topic_table) <- NULL
  out <- list(params = list(), document_topic_table = document_topic_table, 
              log_likelihood = NULL, data = data)
  class(out) <- c("LDA", "list")
  out
}

#' @title Calculate the log likelihood of a VEM LDA model fit
#'
#' @description Imported but updated calculations from topicmodels package, as
#'   applied to Latent Dirichlet Allocation fit with Variational Expectation 
#'   Maximization via \code{\link[topicmodels]{LDA}}. 
#'
#' @details The number of degrees of freedom is 1 (for alpha) plus the number
#'   of entries in the document-topic matrix. The number of observations is 
#'   the number of entries in the document-term matrix.
#'
#' @param object A \code{LDA_VEM}-class object.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @return Log likelihood of the model \code{logLik}, also with \code{df}
#'   (degrees of freedom) and \code{nobs} (number of observations) values.
#'
#' @references 
#'   Buntine, W. 2002. Variational extensions to EM and multinomial PCA. 
#'   \emph{European Conference on Machine Learning, Lecture Notes in Computer 
#'   Science} \strong{2430}:23-34. \href{https://bit.ly/327sltH}{link}.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#'
#'   Hoffman, M. D., D. M. Blei, and F. Bach. 2010. Online learning for 
#'   latent Dirichlet allocation. \emph{Advances in Neural Information 
#'   Processing Systems} \strong{23}:856-864.
#'   \href{https://bit.ly/2LEr5sb}{link}.
#'
#' @examples 
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   r_LDA <- LDA_set(lda_data, topics = 2)   
#'   logLik(r_LDA[[1]])
#'
#' @export
#'


#' @rdname LDA
#'   
#' @export
#'
check_LDA_inputs <- function(data, topics, nreps, 
                                 control){
  check_document_term_table(document_term_table = data$document_term_table)
  check_topics(topics = topics)
  check_nreps(nreps = nreps)
  check_control(control = control)
}




#' @title Select the best LDA model(s) for use in time series
#'
#' @description Select the best model(s) of interest from an
#'   \code{LDA_set} object, based on a set of user-provided functions. The
#'   functions default to choosing the model with the lowest AIC value.
#'
#' @param LDA_models An object of class \code{LDA_set} produced by
#'   \code{\link{LDA_set}}.
#'
#' @param control A \code{list} of parameters to control the running and 
#'   selecting of LDA models. Values not input assume default values set 
#'   by \code{\link{LDA_set_control}}. Values for running the LDAs replace 
#'   defaults in (\code{LDAcontol}, see \code{\link[topicmodels]{LDA}} (but if
#'    \code{seed} is given, it will be overwritten; use \code{iseed} instead).
#'
#' @return A reduced version of \code{LDA_models} that only includes the 
#'   selected LDA model(s). The returned object is still an object of
#'   class \code{LDA_set}.
#'
#' @examples
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   r_LDA <- LDA_set(lda_data, topics = 2, nseeds = 2)  
#'   select_LDA(r_LDA)                       
#'
#' @export
#'
select_LDA <- function(LDA_models = NULL, control = list()){
  if("LDA_set" %in% attr(LDA_models, "class") == FALSE){
    stop("LDA_models must be of class LDA_set")
  }
  control <- do.call("LDA_set_control", control)
  measurer <- control$measurer
  selector <- control$selector  
  lda_measured <- vapply(LDA_models, measurer, 0) %>%
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
#'   (\code{LDA_set}) of the models returned by \code{\link{LDA_set}}.
#'
#' @param mods Fitted models returned from \code{\link[topicmodels]{LDA}}.
#'
#' @param mod_topics Vector of \code{integer} values corresponding to the 
#'   number of topics in each model.
#' 
#' @param mod_seeds Vector of \code{integer} values corresponding to the 
#'   seed used for each model.
#'
#' @return \code{lis} (class: \code{LDA_set}) of LDA models (class: 
#'   \code{LDA_VEM}).
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   topics <- 2
#'   nseeds <- 2
#'   control <- LDA_set_control()
#'   mod_topics <- rep(topics, each = length(seq(2, nseeds * 2, 2)))
#'   iseed <- control$iseed
#'   mod_seeds <- rep(seq(iseed, iseed + (nseeds - 1)* 2, 2), length(topics))
#'   nmods <- length(mod_topics)
#'   mods <- vector("list", length = nmods)
#'   for (i in 1:nmods){
#'     LDA_msg(mod_topics[i], mod_seeds[i], control)
#'     control_i <- prep_LDA_control(seed = mod_seeds[i], control = control)
#'     mods[[i]] <- topicmodels::LDA(document_term_table, k = mod_topics[i], 
#'                      control = control_i)
#'   }
#'   package_LDA_set(mods, mod_topics, mod_seeds)
#' }
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
#' @param LDA_topics \code{integer} value corresponding to the number of 
#'   topics in the model.
#' 
#' @param LDA_rep \code{integer} value corresponding to the replicate number 
#'   for the run of the model.
#'
#' @param control Class \code{LDA_controls} list of control parameters to be
#'   used in \code{LDA} (note that "seed" will be overwritten).
#'
#' @examples
#'   LDA_msg(mod_topics = 4, mod_seeds = 2)
#'
#' @export
#'
LDA_msg <- function(topics, rep, quiet = FALSE){
  check_topics(topics)
  check_nreps(rep)
  topic_msg <- paste0("  ", topics, " topics")
  rep_msg <- paste0(", seed ", rep)
  messageq(paste0(topic_msg, rep_msg), quiet)
}

#' @title Create control list for set of LDA models
#'
#' @description This function provides a simple creation and definition of 
#'   the list used to control the set of LDA models, regardless of how they
#'   are implemented. 
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly.
#'
#' @param measurer,selector Function names for use in evaluation of the LDA
#'   models. \code{measurer} is used to create a value for each model
#'   and \code{selector} operates on the values to choose the model(s) to 
#'   pass on. 
#'
#' @param args \code{list} of arguments to \code{LDA_function}. 
#'
#' @param LDA_function \code{function} used for the LDA model. Defaults to
#'   classical LDA as coded via \code{topicmodels}' 
#'   \code{\link[topicmodels]{LDA}}.
#'
#' @param ... Arguments to be passed to \code{LDA_function} as a 
#'   \code{control} input.
#'
#' @return \code{list} for controlling the LDA model fit.
#'
#' @examples
#'   LDA_set_control()
#'
#' @export
#'
LDA_control <- function(LDA_function = topicmodels_LDA, 
                        LDA_args = list(method = "VEM", seeded = TRUE),
                        measurer_function = AIC,
                        measurer_args = list(),
                        selector_function = min,
                        selector_args = list(), 
                        nsubsets = 1,
                        subset_rule = NULL,
                        soften = TRUE, 
                        quiet = FALSE){
  list(LDA_function = LDA_function,  LDA_args = LDA_args, 
       measurer_function = measurer_function, measurer_args = measurer_args, 
       selector_function = selector_function, selector_args = selector_args,
       nsubsets = nsubsets, subset_rule = subset_rule,
       soften = soften, quiet = quiet)
}

