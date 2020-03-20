##
#
#  working in here on documentation and follow through
#  also make sure when no change point that it still says estimating regressor



TS <- function(LDAs, data, formulas = ~ 1, nchangepoints = 0, 
               timename = "time", weights = NULL, control = list()){
  control <- do.call("TS_control", control)
  messageq("----- Time Series Analyses -----", control$quiet)
  TSs <- prep_TS_models(LDAs = LDAs, data = data, formulas = formulas,
                        nchangepoints = nchangepoints, timename = timename,
                        weights = weights, control = control)
  nTS <- length(TSs)
  for (i in 1:nTS){
    TSs[[i]] <- sequential_TS(TS = TSs[[i]], control = control)
  }
  selected_TSs <- select_TS(TSs = TSs, control = control)
  package_TS(selected_TSs = selected_TSs, TSs = TSs, control = control)
}



prep_TS_models <- function(LDAs, data, formulas = ~ 1, nchangepoints = 0, 
                           timename = "time", weights = NULL, 
                           control = list()){

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
    ts_data$train$ts_data$gamma <- lda$document_topic_matrix

    ts_data$test$ts_data <- lda$data$test$document_covariate_table
    ts_data$test$ts_data$gamma <- lda$test_document_topic_matrix

    weights <- iftrue(weights, 
                      document_weights(lda$data$train$document_term_table))

    TSs[[i]] <- list(data = ts_data,
                     data_subset = lda[["data_subset"]],
                     formula = tab$formula[[i]],
                     nchangepoints = tab$nchangepoints[i], 
                     weights = weights,
                     timename = timename,
                     response = control$response,
                     topics = lda$topics, rep = lda$rep)
  }
  name_tab <- data.frame(paste("LDA", tab[ , 1]), 
                         paste(",", tab[ , 2]),
                         paste(",", tab[ , 3], "change points"))
  names(TSs) <- apply(name_tab, 1, paste0, collapse = "")
  TSs

}



select_TS <- function(TSs, control = list()){

  vals <- measure_TS(TSs = TSs, control = control)
  fun <- control$selector
  args <- update_list(control$selector_args, x = vals)
  selection <- do.call(what = fun, args = args)
  TSs[selection]  
}

measure_TS <- function(TSs, control = list()){
  fun <- control$measurer
  args <- control$measurer_args
  nTSs <- length(TSs)
  vals <- rep(NA, nTSs)
  for(i in 1:nTSs){
    args <- update_list(args, object = TSs[[i]])
    vals_i <- do.call(what = fun, args = args)
    if(length(vals_i) != 0){
      vals[i] <- vals_i
    }
  }
  vals
}

package_TS <- function(selected_TSs, TSs, control = list()){
  out <- list(selected_TSs = selected_TSs, TSs = TSs, control = control)
  class(out) <- c("TS_set", "list")
  out
}



#' @title Create the controls list for the Time Series model
#'
#' @description This function provides a simple creation and definition of a
#'   list used to control the time series model fit occurring within 
#'   \code{\link{TS}}. 
#'
#' @param response \code{character} element indicating the response variable 
#'   used in the time series. \cr \cr
#'   Must have a corresponding \code{<response>_TS} function.
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
#' @return \code{list}, with named elements corresponding to the arguments.
#'
#' @examples
#'   TS_control()
#'
#' @export
#'
TS_control <- function(model = sequential_TS,
                       model_args = list(),
                       response = "multinom",
                       response_args = list(),
                       method = "ldats_classic",
                       method_args = ldats_classic_control(),
                       summary_prob = 0.95,
                       measurer = AIC,
                       measurer_args = list(),
                       selector = which.min,
                       selector_args = list(), 
                       soften = TRUE, 
                       quiet = FALSE, ...){
  list(response = response, method = method, method_args = method_args, 
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
#'   the time series) and the range of the time series. If there are no 
#'   change points (i.e. \code{change points} is \code{NULL}, there is still a
#'   single chunk defined by the start and end of the time series.
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables, but specifically here containing the column indicated
#'   by the \code{timename} input. 
#'
#' @param change points Numeric vector indicating locations of the change 
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
#' @param change points Numeric vector indicating locations of the change 
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



