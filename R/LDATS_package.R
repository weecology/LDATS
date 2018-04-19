#' @importFrom magrittr %>%
#' @importFrom foreach %dopar%
#' @importFrom grDevices devAskNewPage rgb
#' @importFrom graphics axis mtext par plot points rect text
#' @importFrom methods is
#' @importFrom stats as.formula logLik rgeom runif
#' @importFrom utils globalVariables
#'
globalVariables("i")
globalVariables("model")
globalVariables("x")

#' @title Performs two-stage LDA-timeseries analyses
#'
#' @description This package is designed to analyze multivariate time series
#'   using a combination of Latent Dirichlet Allocation and multinomial 
#'   regression.
#'
#' @name LDATS
#' @docType package
#' @keywords package
#'
NULL