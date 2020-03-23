
#' @title Conduct a Bayesian compositional Time Series analysis 
#'
#' @description Analyze compositional Time Series models using Bayesian
#'   sampling methods. \cr \cr
#'   \code{TS} is the main interface function for the LDATS application
#'     of Bayesian change point Time Series analyses (Christensen 
#'     \emph{et al.} 2018). \cr \cr
#'   \code{prepare_TS} pre-prepares the TS model objects for simpler 
#'     use within the subfunctions. \cr \cr
#'   \code{check_TS} ensures that the inputs are proper. 
#'     See \code{\link{check_LDAs}}, 
#'     \code{\link{check_document_covariate_table}}, 
#'     code{\link{check_formulas}}, \code{\link{check_nchangepoints}}, 
#'     \code{\link{check_timename}}, \code{\link{check_weights}}, 
#'     and \code{\link{check_control}} for specifics. \cr \cr
#'   \code{TS_control} defines and creates the control \code{list} for the TS
#'     model running. \cr \cr
#'   \code{run_TS} runs (via \code{\link{TS_call}}) all TS models
#'     as set up by \code{prep_TS_models}. \cr \cr
#'   \code{TS_call} runs (via \code{\link{do.call}}) a single TS model
#'     as set up by \code{prep_TS_models}. \cr \cr
#'   \code{TS_msg} produces a model-running message if desired. \cr \cr
#'   \code{measure_TS} determines the fit value used to select among the 
#'     models. \cr \cr
#'   \code{select_TS} chooses the best model(s) of interest based on their
#'     measured values and the selector function. \cr \cr
#'   \code{package_TS} sets the class and names the elements of the results
#'     \code{list} from \code{\link{TS_call}} applied to the 
#'     combination of TS models requested for the LDA model(s) input.
#'
#' @details For a (potentially subset) dataset consisting of proportions of
#'   topics across multiple documents in a corpus 
#'   \enumerate{
#'     \item Conduct multiple compositional Bayesian TS models 
#'       (e.g., changepoint softmax regression; Ripley 1996, Venables 
#'       and Ripley 2002, Western and Kleykamp 2004, Bishop 2006, Ruggieri 
#'       2013) via a generalized linear modeling approach (McCullagh and 
#'       Nelder 1989) and using parallel tempering Markov Chain Monte Carlo
#'      (ptMCMC) methods (Earl and Deem 2005),
#'     \item Select from the TS model results to pick those used to summarize
#'       the whole model, and
#'     \item Package the results.
#'   }
#'
#' @param formulas Vector of \code{\link[stats]{formula}}(s) defining the 
#'   regression between the change points. Any predictor variable included 
#'   must also be a column in \code{data} and any (compositional) response 
#'   variable must be a set of columns in \code{data}. \cr
#'   Each element (formula) in the vector is evaluated for each number of 
#'   change points and each LDA model.
#'
#' @param nchangepoints \code{integer}-conformable vector corresponding to the 
#'   number of change points to include in the models. 0 is valid (corresponds
#'   to no change points, so a singular time series model) and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segmentation of the 
#'   time series into chunks fit with separate models dictated by 
#'   \code{formula}. \cr
#'   Each element in the vector is the number of change points 
#'   used to segment the data for each formula (entry in \code{formulas}) 
#'   component of the TS model, for each selected LDA model.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior. 
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{\link{TS_call}} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of, e.g., \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using \code{document_weights}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'  
#' @param TS,TSs time series model \code{list} (\code{TS}) or a \code{list} 
#'   of many time series model \code{list}s (\code{TSs}).
#'
#' @param selected_TSs \code{list} of selected time series model \code{list}s.
#'
#' @param response \code{character} element indicating the response variable 
#'   used in the time series. \cr \cr
#'   Must have a corresponding \code{<response>_TS} function.
#'
#' @param response_args \code{list} of (named) arguments to be used in 
#'   \code{response} via \code{\link{do.call}}. 
#'   \cr \cr
#'   Could be managed via a \code{<reponse>_TS_control} function like
#'   \code{\link{multinom_TS_control}}.
#'
#' @param summary_prob Probability used for summarizing the posterior 
#'   distributions (via the highest posterior density interval, see
#'   \code{\link[coda]{HPDinterval}}).
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param soften \code{logical} indicator of whether the model should error 
#'   softly or if errors should trigger a full-stop to the pipeline.
#'
#' @param method \code{function} used to drive the sampler of the TS
#'   models; \code{method} defines and operates the computational procedure.
#'   \cr \cr
#'   Current pre-built options include \code{\link{ldats_classic}}.
#'
#' @param method_args \code{list} of (named) arguments to be used in 
#'   \code{method} via \code{\link{do.call}}. 
#'   \cr \cr
#'   Could be managed via a \code{<method>_control} function like
#'   \code{\link{ldats_classic_control}}.
#'
#' @param measurer \code{function} used in evaluation of the TS
#'   models; \code{measurer} creates a value for each model.
#'
#' @param measurer_args \code{list} of (named) arguments to be used in 
#'   \code{measurer} via \code{\link{do.call}}. 
#'
#' @param selector \code{function} usde in evaluation of the TS
#'   models; \code{selector} operates on the values to choose the models. 
#'
#' @param selector_args \code{list} of (named) arguments to be used in 
#'   \code{selector} via \code{\link{do.call}}. 
#'
#' @param ... Not passed along to the output, rather included to allow for
#'   automated removal of unneeded controls.
#'
#' @return 
#'   \code{TS},\code{pacakage_TS}: class \code{TS_set} \code{list} of both 
#'     selected and all results from \code{\link{TS_call}} applied for 
#'     each model on each LDA model input as well as the control \code{list} 
#'     used to fit the model. \cr \cr
#'   \code{prepare_TS}: \code{list} of \code{list}s, each of which is a
#'     preliminary model object for a Time Series model fit. \cr \cr
#'   \code{check_TS}: an error message is thrown if any input is improper, 
#'     otherwise \code{NULL}.
#'   \code{TS_control}: \code{list} of named control elements for
#'     model fitting.
#'   \code{measure_TS}: \code{vector} of values corresponding to the model
#'     evaluations. \cr \cr
#'   \code{select_TS}: \code{list} of selected models' \code{list}s. \cr \cr
#'   \code{run_TS}: \code{TS_set} \code{list} of model results from all
#'     runs of a \code{<model>} function, such as 
#'     \code{\link{topicmodels_TS}}. \cr \cr
#'   \code{TS_call}: \code{TS} \code{list} of model results from a single
#'     run of a \code{<model>} function, such as 
#'     \code{\link{sequential_TS}}. \cr \cr
#'   \code{TS_msg}: a message is produced.
#'
#' @references 
#'   Bishop, C. M. 2006. \emph{Pattern Recognition and Machine Learning}. 
#'    Springer, New York, NY, USA.
#'
#'   Christensen, E., D. J. Harris, and S. K. M. Ernest. 2018.
#'   Long-term community change through multiple rapid transitions in a 
#'   desert rodent community. \emph{Ecology} \strong{99}:1523-1529. 
#'   \href{https://doi.org/10.1002/ecy.2373}{link}.
#'
#'   Earl, D. J. and M. W. Deem. 2005. Parallel tempering: theory, 
#'   applications, and new perspectives. \emph{Physical Chemistry Chemical 
#'   Physics} \strong{7}: 3910-3916.
#'   \href{https://doi.org/10.1039/B509983H}{link}.
#'
#'   McCullagh, P. and J. A. Nelder. 1989. \emph{Generalized Linear Models}.
#'   2nd Edition. Chapman and Hall, New York, NY, USA.
#'
#'   Ripley, B. D. 1996. \emph{Pattern Recognition and Neural Networks}. 
#'   Cambridge University Press, Cambridge, UK.
#'
#'   Ruggieri, E. 2013. A Bayesian approach to detecting change points in 
#'   climactic records. \emph{International Journal of Climatology}
#'   \strong{33}:520-528.
#'   \href{https://doi.org/10.1002/joc.3447}{link}.
#'
#'   Venables, W. N. and B. D. Ripley. 2002. \emph{Modern and Applied
#'   Statistics with S}. Fourth Edition. Springer, New York, NY, USA.
#'
#'   Western, B. and M. Kleykamp. 2004. A Bayesian change point model for 
#'   historical time series analysis. \emph{Political Analysis}
#'   \strong{12}:354-374.
#'   \href{https://doi.org/10.1093/pan/mph023}{link}.
#'
#' @name TS
#'


#' @rdname TS
#'
#' @export
#'
TS <- function(LDAs, formulas = ~ 1, nchangepoints = 0, 
               timename = "time", weights = NULL, control = list()){
  TSs <- prepare_TS(LDAs = LDAs, formulas = formulas,
                    nchangepoints = nchangepoints, timename = timename,
                    weights = weights, control = control)
  TSs <- run_TS(TSs = TSs)
  TSs <- package_TS(TSs = TSs)
  TSs
}

#' @rdname TS
#'
#' @export
#'
check_TS <- function(LDAs, formulas = ~ 1, nchangepoints = 0, 
               timename = "time", weights = NULL, control = list()){
  check_LDAs(LDAs = LDAs)
  check_document_covariate_table(LDAs = LDAs)
  check_formulas(LDAs = LDAs, formulas = formulas)
  check_nchangepoints(nchangepoints = nchangepoints)
  check_timename(LDAs = LDAs, timename = timename)
  check_weights(weights = weights)
  check_control(control = control)
}

#' @rdname TS
#'
#' @export
#'
run_TS <- function(TSs){
  nTS <- length(TSs)
  for (i in 1:nTS){
    TSs[[i]] <- TS_call(TS = TSs[[i]])
  }
  TSs
}


#' @rdname TS
#'
#' @export
#'
TS_call <- function(TS){
  TS_msg(TS = TS)
  fun <- TS$control$model
  args <- update_list(TS$control$model_args, TS = TS)
  soft_call(what = fun, args = args, soften = TS$control$soften)
}



#' @rdname TS
#'
#' @export
#'
TS_msg <- function(TS){
  subset_msg <- paste0("  - data subset ", TS$data_subset)
  topic_msg <- paste0(", ", TS$topics, " topics")
  rep_msg <- paste0(", replicate ", TS$rep)

  formula_msg <- paste0(", ", deparse(TS$formula))
  nchangepoints <- TS$nchangepoints
  txt <- ifelse(nchangepoints == 1, " change point", " change points")
  changepoints_msg <- paste0(", ", nchangepoints, txt)
  msg <- paste0(subset_msg, topic_msg, rep_msg, formula_msg, changepoints_msg)
  messageq(msg, TS$control$quiet)
}


#' @rdname TS
#'
#' @export
#'
prepare_TS <- function(LDAs, formulas = ~ 1, nchangepoints = 0, 
                           timename = "time", weights = NULL, 
                           control = list()){
  check_TS(LDAs = LDAs, formulas = formulas, nchangepoints = nchangepoints , 
           timename = timename, weights = weights, control = control)
  control <- do.call("TS_control", control)
  messageq("----- Time Series Analyses -----", control$quiet)
  if (!is(formulas, "list")) {
    if (is(formulas, "formula")) {
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain formula(s)")
    }
  } else if (!all(vapply(formulas, is, TRUE, "formula"))) {
      stop("formulas does not contain all formula(s)")
  }

  formulas2 <- formulas
  for (i in seq_along(formulas)) {
    tformula <- paste(as.character(formulas[[i]]), collapse = "")
    formulas2[[i]] <- as.formula(paste("gamma", tformula))
  }
  formulas <- formulas2
  nmods <- length(LDAs[[1]])
  mods <- 1:nmods
  tab <- expand.grid(mods, formulas, nchangepoints, stringsAsFactors = FALSE)
  colnames(tab) <- c("LDA", "formula", "nchangepoints") 
  
  nTSs <- NROW(tab)
  TSs <- vector("list", length = nTSs)
  for(i in 1:nTSs){
    lda <- LDAs[[1]][[tab$LDA[i]]]
 
    ts_data <- lda$data
    ts_data$train$ts_data <- lda$data$train$document_covariate_table
    ts_data$train$ts_data$gamma <- lda$document_topic_table

    ts_data$test$ts_data <- lda$data$test$document_covariate_table
    ts_data$test$ts_data$gamma <- lda$test_document_topic_table

    weights <- iftrue(weights, 
                      document_weights(lda$data$train$document_term_table))

    TSs[[i]] <- list(data = ts_data,
                     data_subset = lda[["data_subset"]],
                     formula = tab$formula[[i]],
                     nchangepoints = tab$nchangepoints[i], 
                     weights = weights,
                     timename = timename,
                     control = control,
                     topics = lda$topics, replicate = lda$replicate)
  }
  name_tab <- data.frame(paste("LDA", tab[ , 1]), 
                         paste(",", tab[ , 2]),
                         paste(",", tab[ , 3], "change points"))
  names(TSs) <- apply(name_tab, 1, paste0, collapse = "")
  TSs

}



#' @rdname TS
#'
#' @export
#'
package_TS <- function(TSs){
  selected_TSs <- select_TS(TSs = TSs)
  out <- list(selected_TSs = selected_TSs, TSs = TSs)
  class(out) <- c("TS_set", "list")
  out
}

#' @rdname TS
#'
#' @export
#'
select_TS <- function(TSs){

  vals <- measure_TS(TSs = TSs)
  fun <- TSs[[1]]$control$selector
  args <- update_list(TSs[[1]]$control$selector_args, x = vals)
  args[names(args) == ""] <- NULL
  selection <- do.call(what = fun, args = args)
  TSs[selection]  
}

#' @rdname TS
#'
#' @export
#'
measure_TS <- function(TSs){

  nTSs <- length(TSs)
  vals <- rep(NA, nTSs)
  for(i in 1:nTSs){
    fun <- TSs[[i]]$control$measurer
    args <- TSs[[i]]$control$measurer_args
    args <- update_list(args, object = TSs[[i]])
    args[names(args) == ""] <- NULL
    vals_i <- do.call(what = fun, args = args)
    if(length(vals_i) != 0){
      vals[i] <- vals_i
    }
  }
  vals
}


#' @rdname TS
#'
#' @export
#'
TS_control <- function(model = sequential_TS,
                       model_args = list(control = sequential_TS_control()),
                       response = multinom_TS,
                       response_args = list(control = multinom_TS_control()),
                       method = ldats_classic,
                       method_args = list(control = ldats_classic_control()),
                       summary_prob = 0.95,
                       measurer = AIC,
                       measurer_args = list(NULL),
                       selector = which.min,
                       selector_args = list(NULL), 
                       soften = TRUE, 
                       quiet = FALSE, ...){
  list(model = model, model_args = model_args,
       response = response, response_args = response_args,
       method = method, method_args = method_args, 
       measurer = measurer, measurer_args = measurer_args, 
       selector = selector, selector_args = selector_args,
       summary_prob = summary_prob, soften = soften, quiet = quiet)
}



#' @title Print a Time Series model 
#'
#' @description Convenience function to print only the most important 
#'   components of a \code{TS}-class object fit by 
#'   \code{\link{sequential_TS}}.
#'
#' @param x Class \code{TS_fit} object to be printed.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @return The non-hidden parts of \code{x} are printed and returned 
#'   invisibly as a \code{list}.
#'
#' @export 
#'
print.TS <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !(names(x) %in% hid)
  print(x[notHid])
}


#' @title Determine the log likelihood of a Time Series model
#'
#' @description Convenience function to extract and format the log likelihood
#'   of a \code{TS}-class object fit by \code{\link{sequential_TS}}.
#'
#' @param object Class \code{TS} object to be evaluated.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @return Log likelihood of the model \code{logLik}, also with \code{df}
#'   (degrees of freedom) and \code{nobs} (number of observations) values.
#'
#' @export
#'
logLik.TS <- function(object, ...){
  val <- object$logLik
  attr(val, "df") <- object$nparams
  attr(val, "nobs") <- nrow(object$data)
  class(val) <- "logLik"
  val
}

#' @title Prepare the time chunk table for a change point Time Series model
#'
#' @description Creates the table containing the start and end times for each
#'   chunk within a time series, based on the change points (used to break up
#'   the time series) and the range of the time series. \cr \cr
#'   If there are no change points (i.e. \code{changepoints = NULL}, there is
#'   still a single chunk defined by the start and end of the time series.
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables, but specifically here containing the column indicated
#'   by the \code{timename} input. 
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points. Must be conformable to \code{integer} values. 
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @return \code{data.frame} of \code{start} and \code{end} times (columns)
#'   for each chunk (rows).
#'
#' @export 
#'
prep_chunks <- function(data, changepoints = NULL, timename = "time"){
  start <- c(min(data[ , timename]), changepoints + 1)   
  end <- c(changepoints, max(data[ , timename])) 
  data.frame(start, end)
}


#' @title Verify the change points of a time series model
#'
#' @description Verify that a time series can be broken into a set 
#'   of chunks based on input change points. 
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables.
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points. Must be conformable to \code{integer} values. 
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @return Logical indicator of the check passing \code{TRUE} or failing
#'   \code{FALSE}.
#'
#' @export 
#'
verify_changepoint_locations <- function(data, changepoints = NULL, 
                                     timename = "time"){
  if (is.null(changepoints)){
    return(TRUE)
  }
  first_time <- min(data[ , timename])
  last_time <- max(data[ , timename])
  time_check <- any(changepoints <= first_time | changepoints >= last_time)
  sort_check <- is.unsorted(changepoints, strictly = TRUE)
  !(time_check | sort_check)
}



