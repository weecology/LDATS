#' @title Run a full set of Latent Dirichlet Allocations and Time 
#'   Series models
#'
#' @description Conduct a complete LDATS analysis (Christensen 
#'   \emph{et al.} 2018), including running a suite of Latent Dirichlet
#'   Allocation (LDA) models (Blei \emph{et al.} 2003, Grun and Hornik 2011) 
#'   via \code{\link{LDA_set}}, selecting LDA model(s) via 
#'   \code{\link{select_LDA}}, running a complete set of Bayesian Time Series
#'   (TS) models (Western and Kleykamp 2004) via \code{\link{TS_on_LDA}} on
#'   the chosen LDA model(s), and selecting the best TS model via 
#'   \code{\link{select_TS}}. \cr \cr
#'   \code{conform_LDA_TS_data} converts the \code{data} input to
#'   match internal and sub-function specifications. \cr \cr
#'   \code{check_LDA_TS_inputs} checks that the inputs to 
#'   \code{LDA_TS} are of proper classes for a full analysis.
#' 
#' @param data Either a document term table or a list including at least
#'   a document term table (with the word "term" in the name of the element)
#'   and optionally also a document covariate table  (with the word 
#'   "covariate" in the name of the element). 
#'   \cr \cr 
#'   The document term table is a table of observation count data (rows: 
#'   documents, columns: terms) that may be a \code{matrix} or 
#'   \code{data.frame}, but must be conformable to a matrix of integers,
#'   as verified by \code{\link{check_document_term_table}}.   
#'   \cr \cr
#'   The document covariate table is a table of associated data (rows: 
#'   documents, columns: time index and covariate options) that may be a
#'   \code{matrix} or \code{data.frame}, but must be a conformable to a data 
#'   table, as verified by \code{\link{check_document_covariate_table}}. Every 
#'   model needs a covariate to describe the time value for each document 
#'   (in whatever units and whose name in the table is input in 
#'   \code{timename}) that dictates the application of the change points. 
#'   \strong{\emph{If a covariate table is not provided, the model assumes the 
#'   observations were equi-spaced in time}}. All covariates named within 
#'   specific models in \code{formulas} must be included. 
#'
#' @param topics Vector of the number of topics to evaluate for each model.
#'   Must be conformable to \code{integer} values.
#'
#' @param nseeds \code{integer} number of seeds (replicate starts) to use for 
#'   each value of \code{topics} in the LDAs. Must be conformable to 
#'   \code{integer} value.
#'
#' @param formulas Vector of \code{\link[stats]{formula}}(s) for the 
#'   continuous (non-change point) component of the time series models. Any 
#'   predictor variable included in a formula must also be a column in the
#'   \code{document_covariate_table}. Each element (formula) in the vector
#'   is evaluated for each number of change points and each LDA model.
#'
#' @param nchangepoints Vector of \code{integer}s corresponding to the number 
#'   of change points to include in the time series models. 0 is a valid input 
#'   corresponding to no change points (\emph{i.e.}, a singular time series
#'   model), and the current implementation can reasonably include up to 6 
#'   change points. Each element in the vector is the number of change points 
#'   used to segment the data for each formula (entry in \code{formulas}) 
#'   component of the TS model, for each selected LDA model.
#' 
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @param weights Optional input for overriding standard weighting for 
#'   documents in the time series. Defaults to \code{TRUE},
#'   translating to an appropriate weighting of the documents
#'   based on the size (number of words) each document (the result of 
#'   \code{\link[topicmodels]{LDA}} is a matrix of proportions, which does not
#'   account for size differences among documents. Alternatively can be 
#'   \code{NULL} for an equal weighting among documents or a \code{numeric}
#'   vector.
#'
#' @param control A \code{list} of parameters to control the running and 
#'   selecting of LDA and TS models. Values not input assume default values 
#'   set by \code{\link{LDA_TS_control}}. 
#' 
#' @param quiet \code{logical} indicator for \code{conform_LDA_TS_data} to 
#'   indicate if messages should be printed.
#'
#' @return 
#'   \code{LDA_TS}: a class \code{LDA_TS} list object including all 
#'   fitted LDA and TS models and selected models specifically as elements 
#'   \code{"LDA models"} (from \code{\link{LDA_set}}),
#'   \code{"Selected LDA model"} (from \code{\link{select_LDA}}), 
#'   \code{"TS models"} (from \code{\link{TS_on_LDA}}), and
#'   \code{"Selected TS model"} (from \code{\link{select_TS}}). \cr \cr
#'   \code{conform_LDA_TS_data}: a data \code{list} that is ready for analyses 
#'   using the stage-specific functions. \cr \cr
#'   \code{check_LDA_TS_inputs}: an error message is thrown if any input is 
#'   improper, otherwise \code{NULL}.
#' 
#' @references 
#'   Blei, D. M., A. Y. Ng, and M. I. Jordan. 2003. Latent Dirichlet
#'   Allocation. \emph{Journal of Machine Learning Research} 
#'   \strong{3}:993-1022.
#'   \href{https://jmlr.csail.mit.edu/papers/v3/blei03a.html}{link}.
#'
#'   Christensen, E., D. J. Harris, and S. K. M. Ernest. 2018.
#'   Long-term community change through multiple rapid transitions in a 
#'   desert rodent community. \emph{Ecology} \strong{99}:1523-1529. 
#'   \href{https://pubmed.ncbi.nlm.nih.gov/29718539/}{link}.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#'
#'   Western, B. and M. Kleykamp. 2004. A Bayesian change point model for 
#'   historical time series analysis. \emph{Political Analysis}
#'   \strong{12}:354-374.
#'   \href{https://www.cambridge.org/core/journals/political-analysis/article/abs/bayesian-change-point-model-for-historical-time-series-analysis/F7D2EDBBC211278EC6C6CB43FE170812}{link}.
#'
#' @examples 
#'   data(rodents)
#' \donttest{
#'   mod <- LDA_TS(data = rodents, topics = 2, nseeds = 1, formulas = ~1,
#'                 nchangepoints = 1, timename = "newmoon")
#' }
#'   conform_LDA_TS_data(rodents)
#'   check_LDA_TS_inputs(rodents, timename = "newmoon")
#'
#' @export
#'
LDA_TS <- function(data, topics = 2, nseeds = 1, formulas = ~ 1,
                   nchangepoints = 0, timename = "time", weights = TRUE, 
                   control = list()){
  check_LDA_TS_inputs(data, topics, nseeds, formulas, nchangepoints,  
                      timename, weights, control)
  control <- do.call("LDA_TS_control", control)
  data <- conform_LDA_TS_data(data, control$quiet)
  dtt <- data$document_term_table
  dct <- data$document_covariate_table
  weights <- iftrue(weights, document_weights(dtt))
  messageq("----Latent Dirichlet Allocation----", control$quiet)
  LDAs <- LDA_set(dtt, topics, nseeds, control$LDA_set_control)
  sel_LDA <- select_LDA(LDAs, control$LDA_set_control)
  messageq("----Time Series Models----", control$quiet)
  TSs <- TS_on_LDA(sel_LDA, dct, formulas, nchangepoints, timename, weights, 
                   control$TS_control)
  sel_TSs <- select_TS(TSs, control$TS_control)
  package_LDA_TS(LDAs, sel_LDA, TSs, sel_TSs)
}

#' @rdname LDA_TS
#'
#' @export
#'
conform_LDA_TS_data <- function(data, quiet = FALSE){
  if(inherits(data, "data.frame") | inherits(data, "matrix")){
    msg <- "covariate table not provided, assuming equi-spaced data"
    messageq(msg, quiet)
    nobs <- nrow(data)
    covariate <- data.frame(time = 1:nobs)
    data <- list(document_term_table = data, 
                 document_covariate_table = covariate)
  } else if(inherits(data, "list")){
    which_term <- grep("term", names(data), ignore.case = TRUE)
    which_covariate <- grep("covariate", names(data), ignore.case = TRUE)
    if(length(which_term) != 1){
      stop("one, and only one, element in `data` can include `term`")
    }
    if (length(which_covariate) == 0){
      msg <- "covariate table not provided, assuming equi-spaced data"
      messageq(msg, quiet)
      nobs <- nrow(data[[which_term]])
      covariate <- data.frame(time = 1:nobs)
      data$document_covariate_table <- covariate    
    } else if(length(which_covariate) > 1){
      stop("at most one element in `data` can include `covariate`")
    }
    names(data)[which_term] <- "document_term_table"
    names(data)[which_covariate] <- "document_covariate_table"
  } else{
    stop("data must be a data.frame, matrix, or list")
  }
  data
}

#' @rdname LDA_TS
#' 
#' @export
#'
check_LDA_TS_inputs <- function(data = NULL,
                              topics = 2, nseeds = 1, formulas = ~ 1, 
                              nchangepoints = 0,
                              timename = "time", 
                              weights = TRUE, 
                              control = list()){
  check_control(control)
  control <- do.call("LDA_TS_control", control)
  data <- conform_LDA_TS_data(data)
  weights <- iftrue(weights, document_weights(data$document_term_table))
  check_document_covariate_table(data$document_covariate_table, 
                               document_term_table = data$document_term_table)
  check_timename(data$document_covariate_table, timename)
  check_formulas(formulas, data$document_covariate_table, control$TS_control)  
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_document_term_table(data$document_term_table)
  check_topics(topics)
  check_seeds(nseeds)
}

#' @title Print the selected LDA and TS models of LDA_TS object
#'
#' @description Convenience function to print only the selected elements of a 
#'   \code{LDA_TS}-class object returned by \code{\link{LDA_TS}}
#'
#' @param x Class \code{LDA_TS} object to be printed.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @return The selected models in \code{x} as a two-element \code{list} with
#'   the TS component only returning the non-hidden components.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   mod <- LDA_TS(data = rodents, topics = 2, nseeds = 1, formulas = ~1,
#'                 nchangepoints = 1, timename = "newmoon")
#'   print(mod)
#' }
#'
#' @export
#'
print.LDA_TS <- function(x, ...){
  print(x[["Selected LDA model"]])
  print(x[["Selected TS model"]])
  list(LDA = x[["Selected LDA model"]], TS = x[["Selected TS model"]])
}

#' @title Package the output of LDA_TS
#'
#' @description Combine the objects returned by \code{\link{LDA_set}},
#'   \code{\link{select_LDA}}, \code{\link{TS_on_LDA}}, and
#'   \code{\link{select_TS}}, name them as elements of the list, and
#'   set the class of the list as \code{LDA_TS}, for the return from
#'   \code{\link{LDA_TS}}.
#'
#' @param LDAs List (class: \code{LDA_set}) of LDA models (class: 
#'   \code{LDA}), as returned by \code{\link{LDA_set}}.
#'
#' @param sel_LDA A reduced version of \code{LDAs} that only includes the 
#'   LDA model(s) selected by \code{\link{select_LDA}}. Still should be of
#'   class \code{LDA_set}.
#'
#' @param TSs Class \code{TS_on_LDA} list of results from \code{\link{TS}} 
#'   applied for each model on each LDA model input, as returned by 
#'   \code{\link{TS_on_LDA}}.
#'
#' @param sel_TSs A reduced version of \code{TSs} (of class \code{TS_fit})
#'   that only includes the TS model chosen via \code{\link{select_TS}}. 
#'
#' @return Class \code{LDA_TS}-class object including all fitted models and 
#'   selected models specifically, ready to be returned from 
#'   \code{\link{LDA_TS}}.
#'
#' @examples 
#' \donttest{
#'   data(rodents)
#'   data <- rodents
#'   control <- LDA_TS_control()              
#'   dtt <- data$document_term_table
#'   dct <- data$document_covariate_table
#'   weights <- document_weights(dtt)
#'   LDAs <- LDA_set(dtt, 2, 1, control$LDA_set_control)
#'   sel_LDA <- select_LDA(LDAs, control$LDA_set_control)
#'   TSs <- TS_on_LDA(sel_LDA, dct, ~1, 1, "newmoon", weights,  
#'                    control$TS_control)
#'   sel_TSs <- select_TS(TSs, control$TS_control)
#'   package_LDA_TS(LDAs, sel_LDA, TSs, sel_TSs)
#' }
#'  
#' @export
#'
package_LDA_TS <- function(LDAs, sel_LDA, TSs, sel_TSs){
  if (!("LDA_set" %in% class(LDAs))){
    stop("LDAs not of class LDA_set")
  }
  if (!("LDA_set" %in% class(sel_LDA))){
    stop("sel_LDA not of class LDA_set")
  }
  if (!("TS_on_LDA" %in% class(TSs))){
    stop("TSs not of class TS_on_LDA")
  }
  if (!("TS_fit" %in% class(sel_TSs))){
    stop("sel_TS not of class TS_fit")
  }

  out <- list("LDA models" = LDAs, "Selected LDA model" = sel_LDA,
              "TS models" = TSs, "Selected TS model" = sel_TSs)
  class(out) <- c("LDA_TS", "list")
  out
}

#' @title Create the controls list for the LDATS model
#'
#' @description Create and define a list of control options used to run the
#'   LDATS model, as implemented by \code{\link{LDA_TS}}.
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly.
#'
#' @param measurer_LDA,selector_LDA Function names for use in evaluation of 
#'   the LDA models. \code{measurer_LDA} is used to create a value for each 
#'   model and \code{selector_LDA} operates on the values to choose the model. 
#'
#' @param iseed \code{integer} initial seed for the LDA model set. 
#'
#' @param ... Additional arguments to be passed to 
#'   \code{\link[topicmodels]{LDA}} as a \code{control} input.
#'
#' @param memoise \code{logical} indicator of whether the multinomial 
#'   functions should be memoised (via \code{\link[memoise]{memoise}}). 
#'   Memoisation happens to both \code{\link{multinom_TS}} and 
#'   \code{\link{multinom_TS_chunk}}.
#'
#' @param response \code{character} element indicating the response variable 
#'   used in the time series. Should be set to \code{"gamma"} for LDATS.
#'
#' @param lambda \code{numeric} "weight" decay term used to set the prior
#'   on the regressors within each chunk-level model. Defaults to 0, 
#'   corresponding to a fully vague prior.
#'
#' @param measurer_TS,selector_TS Function names for use in evaluation of the 
#'   TS models. \code{measurer_TS} is used to create a value for each model
#'   and \code{selector_TS} operates on the values to choose the model. 
#'
#' @param ntemps \code{integer} number of temperatures (chains) to use in the 
#'   ptMCMC algorithm.
#'
#' @param penultimate_temp Penultimate temperature in the ptMCMC sequence.
#'
#' @param ultimate_temp Ultimate temperature in the ptMCMC sequence.
#'
#' @param q Exponent controlling the ptMCMC temperature sequence from the 
#'   focal chain (reference with temperature = 1) to the penultimate chain. 0
#'   (default) implies a geometric sequence. 1 implies squaring before 
#'   exponentiating.
#'
#' @param nit \code{integer} number of iterations (steps) used in the ptMCMC
#'   algorithm.
#'
#' @param magnitude Average magnitude (defining a geometric distribution)
#'   for the proposed step size in the ptMCMC algorithm.
#'
#' @param burnin \code{integer} number of iterations to remove from the 
#'   beginning of the ptMCMC algorithm.
#'
#' @param thin_frac Fraction of iterations to retain, from the ptMCMC. Must be
#'   \eqn{(0, 1]}, and the default value of 1 represents no thinning.
#'
#' @param summary_prob Probability used for summarizing the posterior 
#'   distributions (via the highest posterior density interval, see
#'   \code{\link[coda]{HPDinterval}}) of the TS model.
#'
#' @param seed Input to \code{set.seed} in the time series model for 
#'   replication purposes.
#'
#' @return \code{list} of control \code{lists}, with named elements 
#'   \code{LDAcontrol}, \code{TScontrol}, and \code{quiet}.
#'
#' @examples
#'   LDA_TS_control()
#'
#' @export
#'
LDA_TS_control <- function(quiet = FALSE, measurer_LDA = AIC, 
                           selector_LDA = min, iseed = 2,
                           memoise = TRUE, response = "gamma", lambda = 0, 
                           measurer_TS = AIC, selector_TS = min, ntemps = 6, 
                           penultimate_temp = 2^6, ultimate_temp = 1e10, 
                           q = 0, nit = 1e4, magnitude = 12, burnin = 0, 
                           thin_frac = 1, summary_prob = 0.95, 
                           seed = NULL, ...){

  LDAcontrol <- LDA_set_control(quiet = quiet, measurer = measurer_LDA, 
                                selector = selector_LDA, iseed = iseed, ...)
  TScontrol <- TS_control(memoise = memoise, response = response, 
                          lambda = lambda, measurer = measurer_TS, 
                          selector = selector_TS, ntemps = ntemps, 
                          penultimate_temp = penultimate_temp, 
                          ultimate_temp = ultimate_temp, q = q, 
                          nit = nit, magnitude = magnitude, quiet = quiet, 
                          burnin = burnin, thin_frac = thin_frac, 
                          summary_prob = summary_prob, seed = seed)
  list(LDA_set_control = LDAcontrol, TS_control = TScontrol, quiet = quiet)
}
