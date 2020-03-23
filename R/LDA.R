
#' @title Run a set of Linguistic Decomposition Analysis models
#' 
#' @description Conduct Linguistic Decomposition Analyses. \cr \cr
#'   \code{LDA} provides the main interface for Linguistic Decomposition 
#'     Analysis conducted within the LDATS application of (Christensen 
#'     \emph{et al.} 2018). \cr \cr
#'   \code{prepare_LDA} pre-prepares the LDA model objects for simpler 
#'     use within the subfunctions. \cr \cr 
#'   \code{check_LDA} ensures that the inputs are proper. 
#'     See \code{\link{check_topics}}, \code{\link{check_replicates}},
#'     and \code{\link{check_control}} for specifics. \cr \cr
#'   \code{LDA_control} defines and creates the control list used to fit 
#'     the LDA model. \cr \cr 
#'   \code{run_LDA} runs (via \code{\link{LDA_call}}) all LDA models
#'     as set up by \code{prep_LDA_models}. \cr \cr
#'   \code{LDA_call} runs (via \code{\link{do.call}}) a single LDA model
#'     as set up by \code{prep_LDA_models}. \cr \cr
#'   \code{LDA_msg} produces a model-running message if desired. \cr \cr
#'   \code{measure_LDA} determines the fit value used to select among the 
#'     models. \cr \cr
#'   \code{select_LDA} chooses the best model(s) of interest based on their
#'     measured values and the selector function. \cr \cr
#'   \code{package_LDA} sets the class and names the elements of the results
#'     \code{list} from \code{\link{LDA_call}} applied to the 
#'     combination of TS models requested for the data input(s).
#'
#' @details For a (potentially subset) dataset consisting of counts of words 
#'   across multiple documents in a corpus, 
#'   \enumerate{
#'     \item Conduct multiple Linguistic Decomposition Analysis (LDA) models 
#'       (e.g., Latent Dirichlet Allocation using the Variational Expectation
#'       Maximization (VEM) algorithm; Blei \emph{et al.} 2003, Grun and
#'       Hornik 2011),
#'     \item Select from the LDA model results to pick those used in the Time
#'       Series (TS) models, and
#'     \item Package the results.
#'   }
#'
#' @param data Any of the data structures allowable for LDATS analyses:
#'   \code{matrix} or \code{data.frame} document term table, 
#'   \code{list} of document term and covariate tables, a \code{list} of 
#'   training and test sets of the two tables, or a \code{list} of multiple 
#'   replicate splits of training and test sets of the two tables. \cr
#'   See \code{\link{conform_data}}, which is used to ensure data structure
#'   validity for the desired model.
#'  
#' @param topics Vector of the number of topics to evaluate for each model.
#'   Must be conformable to \code{integer} values.
#'
#' @param replicates Number of replicate starts to use for each 
#'   value of \code{topics}. Must be conformable to \code{integer} value.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   LDA model. Values not input assume defaults set by 
#'   \code{\link{LDA_control}}.
#'
#' @param LDA,LDAs model \code{list} (\code{LDA}) or a \code{list} of LDA 
#'   model \code{list}s (\code{LDAs}).
#'
#' @param model Main LDA \code{function}.
#'
#' @param model_args \code{list} of (named) arguments to be used in 
#'   \code{model} via \code{\link{LDA_call}}. 
#'
#' @param nsubsets Number of data subsets.
#'
#' @param subset_rule \code{function} used to subset the data.
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param soften \code{logical} indicator of whether the model should error 
#'   softly or if errors should trigger a full-stop to the pipeline.
#'
#' @param measurer \code{function} used in evaluation of the LDA
#'   models; \code{measurer} creates a value for each model.
#'
#' @param measurer_args \code{list} of (named) arguments to be used in 
#'   \code{measurer} via \code{\link{do.call}}. 
#'
#' @param selector \code{function} usde in evaluation of the LDA
#'   models; \code{selector} operates on the values to choose the models. 
#'
#' @param selector_args \code{list} of (named) arguments to be used in 
#'   \code{selector} via \code{\link{do.call}}. 
#'
#' @param ... Not passed along to the output, rather included to allow for
#'   automated removal of unneeded controls.
#'
#' @return 
#'   \code{LDA},\code{pacakage_LDA}: class \code{LDA_set} \code{list} of 
#'     both selected and all results from \code{\link{LDA_call}} applied for 
#'     each model on each data input(s) as well as the control \code{list} 
#'     used to fit the model. \cr \cr
#'   \code{prepare_LDA}: \code{list} of \code{list}s, each of which is a
#'     preliminary model object for an LDA model fit. \cr \cr
#'   \code{check_LDA}: an error message is thrown if any input is improper, 
#'     otherwise \code{NULL}.
#'   \code{LDA_control}: \code{list} of controls for the LDA model, with
#'     named elements corresponding to the arguments.
#'   \code{run_LDA}: \code{LDA_set} \code{list} of model results from all
#'     runs of a \code{<model>} function, such as 
#'     \code{\link{topicmodels_LDA}}. \cr \cr
#'   \code{LDA_call}: \code{LDA} \code{list} of model results from a single
#'     run of a \code{<model>} function, such as 
#'     \code{\link{topicmodels_LDA}}. \cr \cr
#'   \code{measure_LDA}: \code{vector} of values corresponding to the model
#'     evaluations. \cr \cr
#'   \code{select_LDA}: \code{list} of selected models' \code{list}s.
#' 
#' @references 
#'   Blei, D. M., A. Y. Ng, and M. I. Jordan. 2003. Latent Dirichlet
#'   Allocation. \emph{Journal of Machine Learning Research} 
#'   \strong{3}:993-1022.
#'   \href{http://jmlr.csail.mit.edu/papers/v3/blei03a.html}{link}.
#'
#'   Christensen, E., D. J. Harris, and S. K. M. Ernest. 2018.
#'   Long-term community change through multiple rapid transitions in a 
#'   desert rodent community. \emph{Ecology} \strong{99}:1523-1529. 
#'   \href{https://doi.org/10.1002/ecy.2373}{link}.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#'
#' @name LDA
#'
#'


#' @rdname LDA
#'
#' @export
#'
LDA <- function(data, topics = 2, replicates = 1, control = list()){
  LDAs <- prepare_LDA(data = data, topics = topics, replicates = replicates, 
                      control = control)
  LDAs <- run_LDA(LDAs = LDAs)
  LDAs <- package_LDA(LDAs = LDAs)
  LDAs
}


#' @rdname LDA
#'
#' @export
#'
check_LDA <- function(topics = 2, replicates = 1, control = list()){
  check_topics(topics = topics)
  check_replicates(replicates = replicates)
  check_control(control = control)
}

#' @rdname LDA
#'
#' @export
#'
prepare_LDA <- function(data, topics = 2, replicates = 1, control = list()){
  check_LDA(topics = topics, replicates = replicates, control = control)
  control <- do.call("LDA_control", control)
  messageq("----- Linguistic Decomposition Analyses -----", control$quiet)
  data <- conform_data(data = data, control = control)
  subsets <- names(data)
  if(length(replicates) < length(topics)){
    replicates <- rep(replicates, length(topics))
  }
  LDA_topics <- rep(topics, replicates)
  LDA_reps <- sequence(replicates)
  LDA_subsets <- rep(subsets, each = length(LDA_reps))  
  LDA_reps <- rep(LDA_reps, length(subsets))
  LDA_topics <- rep(LDA_topics, length(subsets))
  nLDA <- length(LDA_topics)
  LDAs <- vector("list", length = nLDA)
  topicword <- rep(NA, nLDA)
  for(i in 1:nLDA){
    LDAs[[i]] <- list(data = data[[LDA_subsets[i]]], 
                      data_subset = LDA_subsets[i],
                      topics = LDA_topics[i], 
                      replicate = LDA_reps[i],
                      control = control)
    topicword[i] <- ifelse(LDA_topics[i] == 1, "topic", "topics")
  }
  name_tab <- data.frame(paste("data subset", LDA_subsets), 
                         paste(",", LDA_topics, topicword),
                         paste(", replicate", LDA_reps))
  names(LDAs) <- apply(name_tab, 1, paste0, collapse = "")

  LDAs
}



#' @rdname LDA
#'
#' @export
#'
run_LDA <- function(LDAs){
  nLDA <- length(LDAs)
  for (i in 1:nLDA){
    LDAs[[i]] <- LDA_call(LDA = LDAs[[i]])
  }
  LDAs
}


#' @rdname LDA
#'
#' @export
#'
LDA_call <- function(LDA){
  LDA_msg(LDA = LDA)
  fun <- LDA$control$model
  args <- update_list(LDA$control$model_args, LDA = LDA)
  soft_call(what = fun, args = args, soften = LDA$control$soften)
}


#' @rdname LDA
#'
#' @export
#'
LDA_msg <- function(LDA){
  subset_msg <- paste0("  - data subset ", LDA$data_subset)
  topic_word <- ifelse(LDA$topics == 1, " topic", " topics")
  topic_msg <- paste0(", ", LDA$topics, topic_word)
  rep_msg <- paste0(", replicate ", LDA$rep)
  messageq(paste0(subset_msg, topic_msg, rep_msg), LDA$control$quiet)
}


#' @rdname LDA
#'
#' @export
#'
package_LDA <- function(LDAs){
  selected_LDAs <- select_LDA(LDAs = LDAs)
  out <- list(selected_LDAs = selected_LDAs, LDAs = LDAs)
  class(out) <- c("LDA_set", "list")
  out
}


#' @rdname LDA
#'
#' @export
#'
select_LDA <- function(LDAs){
  nLDAs <- length(LDAs)
  maxtopics <- 0
  for(i in 1:nLDAs){
    maxtopics <- max(c(maxtopics, LDAs[[i]]$topics))
  }
  if(maxtopics == 1){
    return(LDAs)
  }
  vals <- measure_LDA(LDAs = LDAs)
  fun <- LDAs[[1]]$control$selector
  args <- update_list(LDAs[[1]]$control$selector_args, x = vals)
  args[names(args) == ""] <- NULL
  selection <- do.call(what = fun, args = args)
  LDAs[selection]  
}

#' @rdname LDA
#'
#' @export
#'
measure_LDA <- function(LDAs){
  nLDAs <- length(LDAs)
  vals <- rep(NA, nLDAs)
  for(i in 1:nLDAs){
    fun <- LDAs[[i]]$control$measurer
    args <- LDAs[[i]]$control$measurer_args
    args <- update_list(args, object = LDAs[[i]])
    args[names(args) == ""] <- NULL
    vals_i <- do.call(what = fun, args = args)
    if(length(vals_i) != 0){
      vals[i] <- vals_i
    }
  }
  vals
}


#' @rdname LDA
#'
#' @export
#'
LDA_control <- function(model = topicmodels_LDA, 
                        model_args = list(method = "VEM", seeded = TRUE),
                        measurer = AIC,
                        measurer_args = list(NULL),
                        selector = which.min,
                        selector_args = list(NULL), 
                        nsubsets = 1,
                        subset_rule = NULL,
                        soften = TRUE, 
                        quiet = FALSE, ...){
  list(model = model,  model_args = model_args, 
       measurer = measurer, measurer_args = measurer_args, 
       selector = selector, selector_args = selector_args,
       nsubsets = nsubsets, subset_rule = subset_rule,
       soften = soften, quiet = quiet)
}



#' @title Determine the AIC of a Linguistic Decomposition Analysis
#'   model
#'
#' @description Convenience function to extract and format the AIC
#'   of a \code{LDA}-class object fit by \code{\link{LDA_call}}.
#'
#' @param object Class \code{LDA} object to be evaluated.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @param k Per-parameter numeric penalty.
#'
#' @return AIC of the model.
#'
#' @export
#'
AIC.LDA <- function(object, ..., k = 2){
  lls <- logLik(object)
  -2 * as.numeric(lls) + k * attr(lls, "df")
}


#' @title Determine the log likelihood of a Linguistic Decomposition Analysis
#'   model
#'
#' @description Convenience function to extract and format the log likelihood
#'   of a \code{LDA}-class object fit by \code{\link{LDA_call}}.
#'
#' @param object Class \code{LDA} object to be evaluated.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @return Log likelihood of the model \code{logLik}, also with \code{df}
#'   (degrees of freedom) and \code{nobs} (number of observations) values.
#'
#' @export
#'
logLik.LDA <- function(object, ...){
  object$log_likelihood
}

