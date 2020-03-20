ilr_TS <- function(data, formula, changepoints = NULL, 
                   timename = "time", weights = NULL, 
                   control = list()){

  if (!verify_changepoint_locations(data, changepoints, timename)){
    out <- list("chunk models" = NA, "logLik" = -Inf, "chunks" = NA)
    class(out) <- c("ilr_TS_fit", "list")
    return(out)
  }

  TS_chunk_memo <- memoise_fun(ilr_TS_chunk, control$memoise)

  chunks <- prep_chunks(data, changepoints, timename)
  nchunks <- nrow(chunks)
  fits <- vector("list", length = nchunks)
  for (i in 1:nchunks){
    fits[[i]] <- TS_chunk_memo(data, formula, chunks[i, ], timename, weights,
                               control)
  }
  package_chunk_fits2(chunks, fits) # temporarily2
}


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

 # temporarily2

package_chunk_fits2 <- function(chunks, fits){
  nchunks <- nrow(chunks)
  chunk_times <- paste0("(", chunks[ , "start"], " - ", chunks[ , "end"], ")")
  names(fits) <- paste("chunk", 1:nchunks, chunk_times, "model")
  ll <- sum(vapply(fits, logLik, 0))
  out <- list("chunk models" = fits, "logLik" = ll, "chunks" = chunks)
  class(out) <- c("ilr_TS_fit", "list")
  out
}


logLik.ilr_TS_fit <- function(object, ...){
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