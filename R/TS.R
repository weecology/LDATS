
TS <- function(gamma, document_covariate_table, timename, formula, 
               nchangepoints, weights, control = TS_controls_list()){

  ts_memo <- memoise_fun(multinom_ts, control$memoise)


}






TS_controls_list <- function(memoise = TRUE){
  out <- list(memoise = memoise)
  class(out) <- c("TS_controls", "list")
  out
}


#' @title Logical control on whether or not to memoise
#'
#' @description This function provides a simple, logical toggle control on
#'   whether the function \code{fun} should be memoised via
#'   \code{\link[memoise]{memoise}} or not. 
#'
#' @param fun Function name to (potentially) be memoised.
#'
#' @param memoise_fun Logical value indicatiing if \code{fun} should be 
#'   memoised.
#'
#' @return function, memoised if desired
#'
#' @export
#'
memoise_fun <- function(fun, memoise_fun){
  if (memoise_fun){
    fun <- memoise(fun)
  }
  fun
}