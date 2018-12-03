#' @title Run a full set of Latent Dirichlet Allocation models then Time 
#'   Series models
#'
#' @description This function runs a complete suite of Latent Dirichlet
#'   Allocation (LDA) models (via \code{\link{LDA_set}}), selects the choice
#'   LDA model(s) (via \code{\link{select_LDA}}), runs a complete set of Time 
#'   Series (TS) models (via \code{\link{TS_on_LDA}}) on the choice LDA
#'   model(s), then selects the best TS model (via \code{select_TS}). 
#' 
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})). May be a class 
#'   \code{matrix} or \code{data.frame} but must be conformable to
#'   a code of integers. This table is a document-level summary of the data 
#'   noted as \eqn{w} (the word-level topic identity) in the math description. 
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
#' @param topics Vector of the number of topics to evaluate in the LDAs.
#'
#' @param nseeds Integer number of seeds (replicate starts) to use for each 
#'   value of \code{topics} in the LDAs.
#'
#' @param formulas Vector of \code{formula}(s) for the continuous change. Any 
#'   predictor variable included in a formula must also be a column in the
#'   \code{document_covariate_table}. Each element (formula) in the vector
#'   is evaluated for each number of change points and each LDA model.
#'
#' @param nchangepoints Vector of integers corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. Each element 
#'   (number of change points) in the vector is used to dictate the
#'   segementation of the data  for each continuous model and each LDA model.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description. Defaults to value calculated by 
#'   \code{\link{document_weights}}, but can be set as \code{NULL}.
#' 
#' @param control Class \code{LDA_TS_controls} list that contains 
#'   \code{LDA_controls}, \code{TS_controls}, and the top-level \code{quiet}.
#'
#' @return Class \code{LDA_TS} object including all fitted models and selected 
#'   models specifically. 
#'
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   lda_data <- rodents$document_term_table
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   
#'   mod <- LDA_TS(document_term_table, document_covariate_table,
#'                 topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
#'                 weights = document_weights(document_term_table), 
#'                 control = LDA_TS_controls_list())
#' }
#'
#' @export
#'
LDA_TS <- function(document_term_table, document_covariate_table,
                   topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 0,
                   weights = document_weights(document_term_table), 
                   control = LDA_TS_controls_list()){
  check_LDA_TS_inputs(document_term_table, document_covariate_table,
                      topics, nseeds, formulas, nchangepoints, weights, 
                      control)
  qprint("Latent Dirichlet Allocation", "----", control$quiet)
  LDAs <- LDA_set(document_term_table, topics, nseeds, control$LDA_control)
  sel_LDA <- select_LDA(LDAs, control$LDA_control)
  qprint("Time Series Models", "----", control$quiet)
  TSs <- TS_on_LDA(sel_LDA, document_covariate_table, formulas, nchangepoints,
                   weights, control$TS_control)
  sel_TSs <- select_TS(TSs, control$TS_control)
  package_LDA_TS(LDAs, sel_LDA, TSs, sel_TSs)
}

#' @title Verify that all inputs to LDA_TS are proper
#'
#' @description Determines if the function has a complete proper set of inputs
#'   for running the code within \code{LDA_TS}/. 
#' 
#' @param document_term_table Table of observation count data (rows: 
#'   documents (\eqn{M}), columns: terms (\eqn{V})). May be a class 
#'   \code{matrix} or \code{data.frame} but must be conformable to
#'   a code of integers. This table is a document-level summary of the data 
#'   noted as \eqn{w} (the word-level topic identity) in the math description. 
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
#' @param topics Vector of the number of topics to evaluate in the LDAs.
#'
#' @param nseeds Integer number of seeds (replicate starts) to use for each 
#'   value of \code{topics} in the LDAs.
#'
#' @param formulas Vector of \code{formula}(s) for the continuous change. Any 
#'   predictor variable included in a formula must also be a column in the
#'   \code{document_covariate_table}. Each element (formula) in the vector
#'   is evaluated for each number of change points and each LDA model.
#'
#' @param nchangepoints Vector of integers corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. Each element 
#'   (number of change points) in the vector is used to dictate the
#'   segementation of the data  for each continuous model and each LDA model.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Corresponds to the vector \strong{\eqn{v}} in the math 
#'   description. Defaults to value calculated by 
#'   \code{\link{document_weights}}, but can be set as \code{NULL}.
#' 
#' @param control Class \code{LDA_TS_controls} list that contains 
#'   \code{LDA_controls}, \code{TS_controls}, and the top-level \code{quiet}.
#'
#' @export
#'
check_LDA_TS_inputs <- function(document_term_table, document_covariate_table,
                              topics = 2, nseeds = 1, formulas = ~ 1, 
                              nchangepoints = 0,
                              weights = document_weights(document_term_table), 
                              control = LDA_TS_controls_list()){
  check_document_covariate_table(document_covariate_table, 
                                 document_term_table = document_term_table)
  check_timename(document_covariate_table, control$TS_control$timename)
  check_formulas(formulas, document_covariate_table, control$TS_control)  
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_control(control, "LDA_TS_controls")
  check_document_term_table(document_term_table)
  check_topics(topics)
  check_seeds(nseeds)
}

#' @title Print the selected LDA and TS models of LDA_TS object
#'
#' @description Convenience function to print only the selected elements of a 
#'   \code{LDA_TS}-class object.
#'
#' @param x Class \code{LDA_TS} object to be printed.
#'
#' @param ... Not used, simply included to maintain method compatability.
#'
#' @export
#'
print.LDA_TS <- function(x, ...){
  print(x[["Selected LDA model"]])
  print(x[["Selected TS model"]])
}

#' @title Package the output of LDA_TS
#'
#' @description Combine the objects, name them as elements of the list, and
#'   set the class of the list, for the return from \code{\link{LDA_TS}}.
#'
#' @param LDAs List (class: \code{LDA_set}) of LDA models (class: 
#'   "\code{LDA}").
#'
#' @param sel_LDA A reduced version of \code{LDAs} that only includes the 
#'   selected LDA model(s). The returned object is still an object of
#'   class \code{LDA_set}.
#'
#' @param TSs Class \code{TS_on_LDA} list of results from \code{\link{TS}} 
#'   applied for each model on each LDA model input.
#'
#' @param sel_TSs A reduced version of \code{TSs} that only includes the 
#'   selected TS model. The returned object is still an object of
#'   class \code{TS_fit}.
#'
#' @return Class \code{LDA_TS} object including all fitted models and selected 
#'   models specifically. 
#'
#' @export
#'
package_LDA_TS <- function(LDAs, sel_LDA, TSs, sel_TSs){
  if (!("LDA_set" %in% class(LDAs))){
    stop("LDAs not of class LDA_set")
  }
  if (!("LDA_set" %in% class(sel_LDA))){
    stop("sel_LDA not of class LDA_set")
  }
  if (!("TS_on_LDA" %in% class(TSs))){
    stop("TSs not of class TS_on_LDA")
  }
  if (!("TS_fit" %in% class(sel_TSs))){
    stop("sel_TS not of class TS_fit")
  }

  out <- list("LDA models" = LDAs, "Selected LDA model" = sel_LDA,
              "TS models" = TSs, "Selected TS model" = sel_TSs)
  class(out) <- c("LDA_TS", "list")
  out
}

#' @title Create the controls list for the LDATS model
#'
#' @description This function provides a simple creation and definition of the
#'   list used to control the LDATS model including the LDA and TS model 
#'   specifically.
#'
#' @param LDA_control Named list of control parameters to be used in 
#'   \code{LDA} (note that "seed" will be overwritten).
#'
#' @param TS_control Class \code{TS_controls} list, holding control parameters
#'   for the Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls, generated by 
#'   \code{\link{TS_controls_list}}.
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly.
#' 
#' @return control Class \code{LDA_TS_controls} list that contains 
#'   \code{LDA_controls}, \code{TS_controls}, and the top-level \code{quiet}.
#'
#' @export
#'
LDA_TS_controls_list <- function(TS_control = TS_controls_list(),
                                 LDA_control = LDA_controls_list(),
                                 quiet = FALSE){
  out <- list(LDA_control = LDA_control, TS_control = TS_control, 
              quiet = quiet)
  class(out) <- c("LDA_TS_controls", "list")
  out
}
