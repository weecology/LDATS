#' @title Fit a simplex-based change point Time Series model 
#'
#' @description 
#'   \code{simplex_TS} fits a set of simplex regression models 
#'     (Aitchison 1986, Aitchison \emph{et al.} 2002) to a time 
#'     series of compositional data divided into multiple segments (a.k.a. 
#'     chunks) based on given locations for a set of change points, using
#'     e.g., the isometric log ratio (ILR) transformation
#'     (Egozcue \emph{et al.} 2003, Pawlowsky-Glahn 2003. \cr \cr
#'   \code{simplex_TS_chunk} fits a simplex regression model using, e.g., the
#'     ILR transformation to a defined chunk of time (a.k.a. segment)
#'     \code{[chunk$start, chunk$end]} within a time series. \cr \cr 
#'   \code{simplex_TS_control} defines and creates the control \code{list} for
#'     fitting.
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
#'   each document. When using \code{simplex_TS} in a standard LDATS 
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
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param transformation Ratio \code{function} to use for the transformation 
#'   to the simplex geometry. Options include \code{\link[compositions]{alr}},
#'   \code{\link[compositions]{clr}}, and \code{\link[compositions]{ilr}}.
#'
#' @param ... Not passed along to the output, rather included to allow for
#'   automated removal of unneeded controls.
#'
#' @return 
#'   \code{simplex_TS}: \code{TS_fit} \code{list} of [1] chunk-level model 
#'     fits (\code{"chunk models"}), [2] the total log likelihood across 
#'     all chunks (\code{"logLik"}), and [3] a \code{data.frame} of chunk 
#'     beginning and ending times (with columns \code{"start"} and 
#'     \code{"end"}). \cr \cr
#'   \code{simplex_TS_chunk}: fitted model object for the chunk, 
#'     of class \code{lm}. \cr \cr
#'   \code{simplex_TS_control}: \code{list}, with named elements 
#'     corresponding to response function controls.
#'
#' @references
#'   Aitchison, J. 1986. \emph{The Statistical Analysis of Compositional
#'   Data}. Monographs on Statistics and Applied Probability. Chapman & Hall
#'   Ltd., London, UK.
#'
#'   Aitchison, J, C. Barcelo-Vidal, J.J. Egozcue, and V. Pawlowsky-Glahn.
#'   2002. A consise guide to the algebraic geometric structure of the 
#'   simplex, the sample space for compositional data analysis, Terra Nostra, 
#'   Schriften der Alfred Wegener-Stiftung, 03/2003.
#'
#'   Egozcue J.J., V. Pawlowsky-Glahn, G. Mateu-Figueras and C. Barcelo-Vidal.
#'   2003. Isometric logratio transformations for compositional data analysis. 
#'   \emph{Mathematical Geology}, \strong{35}:279-300.
#'
#'   Pawlowsky-Glahn, V. 2003. Statistical modelling on coordinates. In: 
#'   Thio-Henestrosa, S. and J. A. Martin-Fernandez, Eds. 
#'   \emph{Proceedings of the 1st International Workshop on Compositional Data
#'    Analysis}, Universitat de Girona, ISBN 84-8458-111-X.
#'   \href{link}{http://ima.udg.es/Activitats/CoDaWork03}.
#'
#' @name simplex_TS
#'


#' @rdname simplex_TS 
#'
#' @export 
#'
simplex_TS <- function(data, formula, changepoints = NULL, 
                   timename = "time", weights = NULL, 
                   control = list()){
  control <- do.call("simplex_TS_control", control)
  if (!verify_changepoint_locations(data, changepoints, timename)){
    out <- list("chunk models" = NA, "logLik" = -Inf, "chunks" = NA)
    return(out)
  }

  chunks <- prep_chunks(data, changepoints, timename)
  nchunks <- nrow(chunks)
  fits <- vector("list", length = nchunks)
  for (i in 1:nchunks){
    fits[[i]] <- simplex_TS_chunk(data = data, formula = formula, 
                              chunk = chunks[i, ], timename = timename,
                              weights = weights, control = control)
  }
  package_chunk_fits(chunks, fits)
}

#' @rdname simplex_TS 
#'
#' @export 
#'
simplex_TS_chunk <- function(data, formula, chunk, timename = "time",
                              weights = NULL, control = list()){
  formula <- as.formula(format(formula))
  time_obs <- data[ , timename] 
  chunk_start <- as.numeric(chunk["start"])
  chunk_end <- as.numeric(chunk["end"])
  in_chunk <- time_obs >= chunk_start & time_obs <= chunk_end
  simplex_data <- data[ , !grepl("gamma", colnames(data))]
  simplex_data$gamma <- do.call(control$transformation, 
                                list(x = acomp(data$gamma)))


  fit <- lm(formula, simplex_data, weights, subset = in_chunk)
  fit$timevals <- time_obs[which(in_chunk)]
  fit 
}


#' @rdname simplex_TS 
#'
#' @export 
#'
simplex_TS_control <- function(transformation = ilr, quiet = FALSE, ...){
  list(transformation = transformation, quiet = quiet)
}


#' @title Fit a multinomial change point Time Series model
#'
#' @description 
#'   \code{multinom_TS} fits a set of multinomial regression models (via
#'     \code{\link[nnet]{multinom}}, Venables and Ripley 2002) to a time 
#'     series of data divided into multiple segments (a.k.a. chunks) based on 
#'     given locations for a set of change points. \cr \cr
#'   \code{multinom_TS_chunk} fits a multinomial regression model (via
#'     \code{\link[nnet]{multinom}}, Ripley 1996, Venables and Ripley 2002)
#'     to a defined chunk of time (a.k.a. segment)
#'     \code{[chunk$start, chunk$end]} within a time series. \cr \cr
#'   \code{multinom_TS_control} defines and creates the control \code{list}
#'     for fitting.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}). \cr \cr 
#'   Note that the response variables should be formatted as a 
#'   \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. \cr \cr
#'   See \code{Examples}.
#'
#' @param formula \code{\link[stats]{formula}} defining the regression between
#'   relationship the change points. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (multinomial) response variable must be a set of
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
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link{topicmodels_LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using 
#'   \code{\link{document_weights}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @param lambda \code{numeric} "weight" decay term used to set the prior
#'   on the regressors within each chunk-level model. Defaults to 0, 
#'   corresponding to a fully vague prior.
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param ... Not passed along to the output, rather included to allow for
#'   automated removal of unneeded controls.
#'
#' @return 
#'   \code{multinom_TS}: Object of class \code{TS_fit}, which is a 
#'     \code{list} of [1] chunk-level model fits (\code{"chunk models"}), 
#'     [2] the total log likelihood combined across all chunks 
#'     (\code{"logLik"}), and [3] a \code{data.frame} of chunk beginning and
#'     ending times (with columns \code{"start"} and \code{"end"}). \cr \cr
#'   \code{multinom_TS_chunk}: fitted model object for the chunk, 
#'     of classes \code{multinom} and \code{nnet}. \cr \cr
#'   \code{multinom_TS_control}: \code{list}, with named elements 
#'     corresponding to response function controls.
#'
#' @references
#'   Ripley, B. D. 1996. \emph{Pattern Recognition and Neural Networks}. 
#'   Cambridge University Press, Cambridge, UK.
#'
#'   Venables, W. N. and B. D. Ripley. 2002. \emph{Modern and Applied
#'   Statistics with S}. Fourth Edition. Springer, New York, NY, USA.
#'
#' @name multinom_TS 
#'


#' @rdname multinom_TS 
#'
#' @export 
#'
multinom_TS <- function(data, formula, changepoints = NULL, 
                        timename = "time", weights = NULL, 
                        control = list()){
  control <- do.call("multinom_TS_control", control)
  if (!verify_changepoint_locations(data, changepoints, timename)){
    out <- list("chunk models" = NA, "logLik" = -Inf, "chunks" = NA)
    class(out) <- c("TS_fit", "list")
    return(out)
  }

  chunks <- prep_chunks(data, changepoints, timename)
  nchunks <- nrow(chunks)
  fits <- vector("list", length = nchunks)
  for (i in 1:nchunks){
    fits[[i]] <- multinom_TS_chunk(data = data, formula = formula, 
                                   chunk = chunks[i, ], timename = timename, 
                                   weights = weights, control = control)
  }
  package_chunk_fits(chunks, fits)
}

#' @rdname multinom_TS 
#'
#' @export 
#'
multinom_TS_chunk <- function(data, formula, chunk, timename = "time",
                              weights = NULL, control = list()){
  formula <- as.formula(format(formula))
  time_obs <- data[ , timename] 
  chunk_start <- as.numeric(chunk["start"])
  chunk_end <- as.numeric(chunk["end"])
  in_chunk <- time_obs >= chunk_start & time_obs <= chunk_end
  fit <- multinom(formula, data, weights, subset = in_chunk, trace = FALSE,
                  decay = control$lambda)
  fit$timevals <- time_obs[which(in_chunk)]
  fit 
}

#' @rdname multinom_TS 
#'
#' @export 
#'
multinom_TS_control <- function(lambda = 0, quiet = FALSE, ...){
  list(lambda = lambda, quiet = quiet)
}



#' @title Package the output of the chunk-level TS models into a TS_fit list
#'    
#' @description Takes the list of fitted chunk-level models returned from
#'   a \code{<response>_TS_chunk} function and packages it as a 
#'   \code{TS_fit} object. This involves naming the model fits based 
#'   on the chunk time windows, combining the log likelihood values across the 
#'   chunks, and setting the class of the output object. 
#'
#' @param chunks Data frame of \code{start} and \code{end} times for each 
#'   chunk (row).
#'
#' @param fits List of chunk-level fits returned by \code{TS_chunk_memo},
#'   the memoised version of \code{\link{multinom_TS_chunk}}.
#'
#' @return Object of class \code{TS_fit}, which is a list of [1]
#'   chunk-level model fits, [2] the total log likelihood combined across 
#'   all chunks, and [3] the chunk time data table.
#'
#' @export 
#' 
#'
package_chunk_fits <- function(chunks, fits){
  nchunks <- nrow(chunks)
  chunk_times <- paste0("(", chunks[ , "start"], " - ", chunks[ , "end"], ")")
  names(fits) <- paste("chunk", 1:nchunks, chunk_times, "model")
  ll <- sum(vapply(fits, logLik, 0))
  out <- list("chunk models" = fits, "logLik" = ll, "chunks" = chunks)
  class(out) <- c("TS_fit", "list")
  out
}

#' @title Log likelihood of a TS model (as a TS_fit-class list)
#' 
#' @description Convenience function to simply extract the \code{logLik}
#'   element (and \code{df} and \code{nobs}) from a \code{TS_fit}
#'   object fit by a \code{<response>_TS} function.
#'
#' @param object A \code{TS_fit}-class object.
#'
#' @param ... Not used, simply included to maintain method compatibility.
#'
#' @return Log likelihood of the model, as class \code{logLik}, with 
#'   attributes \code{df} (degrees of freedom) and \code{nobs} (the number of
#'   weighted observations, accounting for size differences among documents).
#'
#' @export
#'
logLik.TS_fit <- function(object, ...){
  ll <- object$logLik
  df <- NA
  nobs <- NA
  if (object$logLik != -Inf){
    nchunks <- nrow(object$chunks)
    dfperchunk <- length(coef(object$"chunk models"[[1]]))
    df <- nchunks - 1 + dfperchunk * nchunks
    nobs <- 0
    for(i in 1:nchunks){
      nobs <- nobs + sum(object$"chunk models"[[i]]$weights)
    }
  }
  structure(ll, df = df, nobs = nobs, class = "logLik")  
}

