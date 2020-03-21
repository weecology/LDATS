
#' @title Run a set of Linguistic Decomposition Analysis models
#' 
#' @description Conduct Linguistic Decomposition Analyses. \cr \cr
#'   \code{LDA} provides the main interface for Linguistic Decomposition 
#'     Analysis conducted within the LDATS application of (Christensen 
#'     \emph{et al.} 2018). \cr \cr
#'   \code{prep_LDA_models} pre-prepares the LDA model objects for simpler 
#'     use within the subfunctions. \cr \cr 
#'   \code{LDA_call} runs (via \code{\link{do.call}}) a single LDA model
#'     as set up by \code{prep_LDA_models}. \cr \cr
#'   \code{LDA_msg} produces a model-running message if desired.
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
#' @param reps Number of replicate starts to use for each 
#'   value of \code{topics}. Must be conformable to \code{integer} value.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   LDA model. Values not input assume defaults set by 
#'   \code{\link{LDA_control}}.
#'
#' @return 
#'   \code{LDA}: class \code{LDA_set} \code{list} of both selected and all
#'     results from \code{\link{LDA_call}} applied for each model on each
#'     data input(s) as well as the control \code{list} used to fit the 
#'     model. \cr \cr
#'   \code{prep_LDA_models}: \code{list} of \code{list}s, each of which is a
#'     preliminary model object for an LDA model fit. \cr \cr
#'   \code{LDA_call}: \code{LDA} \code{list} of model results from a single
#'     run of a \code{<model>} function, such as 
#'     \code{\link{topicmodels_LDA}}. \cr \cr
#'   \code{LDA_msg}: a message is produced.
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
#' @export
#'
LDA <- function(data, topics = 2, reps = 1, control = list()){
  control <- do.call("LDA_control", control)
  messageq("----- Linguistic Decomposition Analyses -----", control$quiet)
  LDAs <- prep_LDA_models(data = data, topics = topics, reps = reps,
                          control = control)
  nLDA <- length(LDAs)
  for (i in 1:nLDA){
    LDAs[[i]] <- LDA_call(LDA = LDAs[[i]], control = control)
  }
  selected_LDAs <- select_LDA(LDAs = LDAs, control = control)
  package_LDA(selected_LDAs = selected_LDAs, LDAs = LDAs, control = control)
}


#' @rdname LDA
#'
#' @export
#'
prep_LDA_models <- function(data, topics = 2, reps = 1, control = list()){
  data <- conform_data(data = data, control = control)
  subsets <- names(data)
  if(length(reps) < length(topics)){
    reps <- rep(reps, length(topics))
  }
  LDA_topics <- rep(topics, reps)
  LDA_reps <- sequence(reps)
  LDA_subsets <- rep(subsets, each = length(LDA_reps))  
  LDA_reps <- rep(LDA_reps, length(subsets))
  LDA_topics <- rep(LDA_topics, length(subsets))
  nLDA <- length(LDA_topics)
  LDAs <- vector("list", length = nLDA)
  for(i in 1:nLDA){
    LDAs[[i]] <- list(data = data[[LDA_subsets[[i]]]], 
                      data_subset = LDA_subsets[[i]],
                      topics = LDA_topics[[i]], rep = LDA_reps[[i]])
  }
  names(LDAs) <- paste0("model_", 1:nLDA)
  LDAs
}


#' @rdname LDA
#'
#' @export
#'
LDA_call <- function(LDA = NULL, control = list()){
  control <- do.call("LDA_control", control)  
  LDA_msg(LDA = LDA, control = control)
  fun <- control$model
  args <- update_list(control$model_args, LDA = LDA)
  soft_call(fun = fun, args = args, soften = control$soften)
}


#' @rdname LDA
#'
#' @export
#'
LDA_msg <- function(LDA, control = list()){
  subset_msg <- paste0("  - data subset ", LDA$data_subset)
  topic_msg <- paste0(", ", LDA$topics, " topics")
  rep_msg <- paste0(", replicate ", LDA$rep)
  messageq(paste0(subset_msg, topic_msg, rep_msg), control$quiet)
}




#' @title Create the controls list for the Time Series model
#'
#' @description This function provides a simple creation and definition of a
#'   list used to control the LDA model fit occurring within 
#'   \code{\link{LDA}}. 
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
#' @return \code{list}, with named elements corresponding to the arguments.
#'
#' @examples
#'   LDA_control()
#'
#' @export
#'
LDA_control <- function(model = topicmodels_LDA, 
                        model_args = list(method = "VEM", seeded = TRUE),
                        measurer = AIC,
                        measurer_args = list(),
                        selector = which.min,
                        selector_args = list(), 
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


#' @title Latent Dirichlet Allocation Linguistic Decomposition Analysis
#'   as conducted via the topicmodels package
#'
#' @description Fit the standard LDATS LDA model (a true Latent Dirichlet
#'   Allocation) using \code{\link[topicmodels]{topicmodels::LDA}}
#'   (Grun and Hornik 2011). \cr 
#'   Default methodology is the Variational
#'   Expectation Maximization routine (VEM) as described by 
#'   Blei \emph{et al.} (2003) and implemented by Grun and Hornik (2011). \cr
#'   If the model is defined to only fit one topic, \code{\link{identity_LDA}}
#'   is used by default.
#'
#' @param LDA A prepared (via \code{\link{prep_LDA_models}} LDA model
#'   \code{list}.
#'
#' @param ... Additional arguments to be passed to 
#'   \code{\link[topicmodels]{topicmodels::LDA}} as a \code{control} input.
#'
#' @param seeded \code{logical} indicator of if the LDA should be a seeded
#'   replicate. 
#'
#' @param method Fitting routine used in 
#'   \code{\link[topicmodels]{topicmodels::LDA}}. Currenlty, only
#'   \code{"VEM"} and \code{"Gibbs"} are supported.
#'
#' @return \code{LDA} \code{list}.
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
#' @export
#'
topicmodels_LDA <- function(LDA, method = "VEM", seeded = TRUE, ...){
  data <- LDA$data
  topics <- LDA$topics 
  rep <- LDA$rep
  data_subset <- LDA$data_subset
  if(topics == 1){
    identity_LDA(LDA)
  } else{
    fun_control <- list(...)
    if(seeded){
      fun_control <- update_list(fun_control, seed = rep * 2)
    }
    mod <- topicmodels::LDA(x = data$train$document_term_table, k = topics, 
                            method = method, control = fun_control)
    mod_ll <- sum(mod@loglikelihood)
    alpha <- tryCatch(as.integer(mod@control@estimate.alpha), 
                      error = function(x){0})
    df <- alpha + length(mod@beta)
    attr(mod_ll, "df") <- df
    attr(mod_ll, "nobs") <- mod@Dim[1] * mod@Dim[2]
    class(mod_ll) <- "logLik"
    out <- list(params = list(alpha = mod@alpha, beta = mod@beta),
                document_topic_matrix = mod@gamma, 
                test_document_topic_matrix = NULL, #not yet available
                log_likelihood = mod_ll, data = data,
                topics = topics, rep = rep, data_subset = data_subset)
    class(out) <- c("LDA", "list")
    out
  } 
}

#' @title Identity Linguistic Decomposition Analysis
#'
#' @description This function acts as an "identity" model, wherein the 
#'   output is functionally the input. This allows for "single-topic" models
#'   that do not actually decompose the data to be included in the model set.
#'
#' @param LDA A prepared (via \code{\link{prep_LDA_models}} LDA model
#'   \code{list}.
#'
#' @return \code{LDA} \code{list} with most components as placeholders.
#' 
#' @export
#'
identity_LDA <- function(LDA){
  data <- LDA$data
  rep <- LDA$rep
  data_subset <- LDA$data_subset
  document_topic_table <- data$train$document_term_table 
  document_topic_table <- document_topic_table / rowSums(document_topic_table)
  colnames(document_topic_table) <- NULL
  out <- list(params = list(), document_topic_table = document_topic_table, 
              log_likelihood = NULL, data = data,
              topics = 1, rep = rep, data_subset = data_subset,
              test_document_topic_matrix = NULL) #not yet available
  class(out) <- c("LDA", "list")
  out
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



#' @title Measure, select, and package the output of a set of Linguistic
#'   Decomposition Analysis models
#'
#' @description 
#'   \code{measure_LDA} determines the fit value used to select among the 
#'     models. \cr \cr
#'   \code{select_LDA} chooses the best model(s) of interest based on their
#'     measured values and the selector function. \cr \cr
#'   \code{package_LDA} sets the class and names the elements of the results
#'     \code{list} from \code{\link{LDA_call}} applied to the 
#'     combination of TS models requested for the data input(s).
#'
#' @param LDAs \code{list} of LDA model \code{list}s.
#'
#' @param selected_LDAs \code{list} of selected LDA model \code{list}s.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   LDA model. Values not input assume defaults set by 
#'   \code{\link{LDA_control}}.
#'
#' @return 
#'   \code{measure_LDA}: \code{vector} of values corresponding to the model
#'     evaluations.
#'   \code{select_LDA}: \code{list} of selected models' \code{list}s.
#'   \code{pacakage_LDA}: class \code{LDA_set} \code{list} of both selected 
#'     and all results from \code{\link{LDA_call}} applied for each model on
#'     each data input as well as the control \code{list} used to fit 
#'     the model.
#'
#' @export
#'
package_LDA <- function(selected_LDAs, LDAs, control = list()){
  out <- list(selected_LDAs = selected_LDAs, LDAs = LDAs, control = control)
  class(out) <- c("LDA_set", "list")
  out
}


#' @rdname package_LDA
#'
#' @export
#'
select_LDA <- function(LDAs = list(), control = list()){
  nLDAs <- length(LDAs)
  maxtopics <- 0
  for(i in 1:nLDAs){
    maxtopics <- max(c(maxtopics, LDAs[[i]]$topics))
  }
  if(maxtopics == 1){
    return(LDAs)
  }
  vals <- measure_LDA(LDAs = LDAs, control = control)
  fun <- control$selector
  args <- update_list(control$selector_args, x = vals)
  selection <- do.call(what = fun, args = args)
  LDAs[selection]  
}

#' @rdname package_LDA
#'
#' @export
#'
measure_LDA <- function(LDAs = list(), control = list()){
  fun <- control$measurer
  args <- control$measurer_args
  nLDAs <- length(LDAs)
  vals <- rep(NA, nLDAs)
  for(i in 1:nLDAs){
    args <- update_list(args, object = LDAs[[i]])
    vals_i <- do.call(what = fun, args = args)
    if(length(vals_i) != 0){
      vals[i] <- vals_i
    }
  }
  vals
}



