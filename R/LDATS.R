#' @importFrom digest digest
#' @importFrom graphics axis mtext par plot points rect text
#' @importFrom grDevices devAskNewPage rgb
#' @importFrom lubridate is.Date
#' @importFrom magrittr %>%
#' @importFrom memoise memoise
#' @importFrom methods is
#' @importFrom nnet multinom
#' @importFrom progress progress_bar
#' @importFrom stats AIC as.formula terms
#' @importFrom topicmodels LDA logLik
#' @importFrom viridis viridis
#'

#' @title Two-stage LDA-TimeSeries analyses
#'
#' @description Performs two-stage analysis of multivariate temporal data
#'   using a combination of Latent Dirichlet Allocation and multinomial 
#'   time series models.
#'
#' @name LDATS
#'
#' @docType package
#'
#' @keywords package
#'
NULL