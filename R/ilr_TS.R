#' @title Fit a simplex-based change point Time Series model using isometric
#'   log ratio
#'
#' @description 
#'   \code{ilr_TS} fits a set of simplex regression models to a time series
#'   of compositional data divided into multiple segments (a.k.a. chunks)
#'   based on given locations for a set of change points, using
#'   the isometric log ratio (ILR) transformation. \cr \cr
#'   \code{ilr_TS_chunk} fits a simplex regression model using the ILR
#'   transformation to a defined chunk of time (a.k.a. segment)
#'   \code{[chunk$start, chunk$end]} within a time series.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the compositional response variable (indicated 
#'   in \code{formula}). \cr \cr 
#'   Note that the response variables should be formatted as a 
#'   \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. \cr \cr
#'   See \code{Examples}.
#'
#' @param formula \code{\link[stats]{formula}} defining the regression between
#'   relationship the change points. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (compositional) response variable must be a set of
#'   columns in \code{data}.
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points. Must be conformable to \code{integer} values.
#'
#' @param chunk Length-2 vector of times: [1] \code{start}, the start time 
#'   for the chunk and [2] \code{end}, the end time for the chunk.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{ilr_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link{topicmodels_LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using 
#'   \code{\link{document_weights}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{ilr_TS}: Object of class \code{ilr_TS_fit}, 
#'   which is a list of [1]
#'   chunk-level model fits (\code{"chunk models"}), [2] the total log 
#'   likelihood combined across all chunks (\code{"logLik"}), and [3] a 
#'   \code{data.frame} of chunk beginning and ending times (\code{"logLik"}
#'   with columns \code{"start"} and \code{"end"}). \cr \cr
#'   \code{ilr_TS_chunk}: fitted model object for the chunk, 
#'   of class \code{lm}.
#'
#' @export 
#'

ilr_TS <- function(data, formula, changepoints = NULL, 
                   timename = "time", weights = NULL, 
                   control = list()){
  control <- do.call("ilr_TS_control", control)
  if (!verify_changepoint_locations(data, changepoints, timename)){
    out <- list("chunk models" = NA, "logLik" = -Inf, "chunks" = NA)
    return(out)
  }

  chunks <- prep_chunks(data, changepoints, timename)
  nchunks <- nrow(chunks)
  fits <- vector("list", length = nchunks)
  for (i in 1:nchunks){
    fits[[i]] <- ilr_TS_chunk(data = data, formula = formula, 
                              chunk = chunks[i, ], timename = timename,
                              weights = weights, control = control)
  }
  package_chunk_fits(chunks, fits)
}

#' @rdname ilr_TS 
#'
#' @export 
#'
ilr_TS_chunk <- function(data, formula, chunk, timename = "time",
                              weights = NULL, control = list()){
  formula <- as.formula(format(formula))
  time_obs <- data[ , timename] 
  chunk_start <- as.numeric(chunk["start"])
  chunk_end <- as.numeric(chunk["end"])
  in_chunk <- time_obs >= chunk_start & time_obs <= chunk_end
  ilr_data <- data[ , !grepl("gamma", colnames(data))]
  ilr_data$gamma <- ilr(acomp(data$gamma))
  fit <- lm(formula, ilr_data, weights, subset = in_chunk)
  fit$timevals <- time_obs[which(in_chunk)]
  fit 
}



#' @title Create the controls list for the ILR Time Series model
#'
#' @description This function provides a simple creation and definition of a
#'   list used to control the time series model fit occurring within 
#'   \code{\link{ilr_TS}}. 
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param ... Not passed along to the output, rather included to allow for
#'   automated removal of unneeded controls.
#'
#' @return \code{list}, with named elements corresponding to the arguments.
#'
#' @export
#'
ilr_TS_control <- function(quiet = FALSE, ...){
  list(quiet = quiet)
}
