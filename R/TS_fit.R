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
