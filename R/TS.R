#' @title Conduct a set of Time Series analyses on LDA models
#'
#' @description This is a wrapper function that expands the main Time Series
#'   analyses function (\code{MTS}) across the LDA models and the Time Series 
#'   models, with respect to both continuous time formulae and the number
#'   of discrete changepoints. This function allows direct passage of the
#'   control parameters for the parallel tempering MCMC through to the 
#'   main Time Series function, \code{MTS}, via the \code{ptMCMC_controls}
#'   argument  
#'
#' @param lda_models List of LDA models (class \code{LDA_list}) or a singular
#'   LDA model (class \code{LDA}).
#'
#' @param document_covariate_table Document covariate table (rows:
#'   documents (\code{M}), columns: time index and covariate options). 
#'   Every model needs a covariate to describe the time value for each
#'   document (in whatever relevant units), whose name in the table is input
#'   via \code{timename}, that dictates the application of the changepoints. 
#'   In addition, the table needs to include as columns all covariates named 
#'   within the specific models described via the argument \code{formula} 
#'   (if desired). Must be a conformable to a data table. 
#'
#' @param timename Character name of the column in the 
#'   \code{document_covariate_table} that contains the time index to use
#'   for assignment of the changepoints (corresponding to the vector 
#'   \strong{\eqn{t}} in the mathematical description of the model). 
#'
#' @param formula Vector of \code{formula}(s) for the continuous change. Any 
#'   predictor variable included in a formula must also be a column in the
#'   \code{document_covariate_table}. Each element (formula) in the vector
#'   is evaluated for each number of change points and each LDA model.
#'
#' @param changepoints Vector of integers corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. Each element 
#'   (number of change points) in the vector is used to dictate the
#'   segementation of the data  for each continuous model and each LDA model.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description.
#'
#' @param ptMCMC_controls Class \code{ptMCMC_controls} list, holding control 
#'   parameters for the parallel tempering Markov Chain Monte Carlo (ptMCMC)
#'   fitting of the MTS models.
#'
#' @return the set of results from MTS
#'
#' @export
#'
TS_set_on_LDA <- function(lda_models, document_covariate_table, timename,
                          formula = ~ 1, changepoints = 0, weights = NULL, 
                          ptMCMC_controls = ptMCMC_controls_list()){

}


ptMCMC_controls_list <- function(){
  out <- list()
  class(out) <- c("ptMCMC_controls", "list")
  out
}