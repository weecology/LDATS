#' @title Conduct a single multinomial Bayesian Time Series analysis 
#'
#' @description This is the main interface function for the LDATS application
#'   of Bayesian change point Time Series analyses (Christensen \emph{et al.}
#'   2018), which extends the model of Western and Kleykamp (2004;
#'   see also Ruggieri 2013) to multinomial (proportional) response data using
#'   softmax regression (Ripley 1996, Venables and Ripley 2002, Bishop 2006) 
#'   using a generalized linear modeling approach (McCullagh and Nelder 1989).
#'   The models are fit using parallel tempering Markov Chain Monte Carlo
#'   (ptMCMC) methods (Earl and Deem 2005) to locate change points and 
#'   neural networks (Ripley 1996, Venables and Ripley 2002, Bishop 2006) to
#'   estimate regressors.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. See \code{Examples}.
#'
#' @param formula \code{\link[stats]{formula}} defining the regression between
#'   relationship the change points. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (multinomial) response variable must be a set of
#'   columns in \code{data}, as verified by \code{\link{check_formula}}.
#'
#' @param nchangepoints \code{integer} corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segementation of the 
#'   time series into chunks fit with separate models dictated by 
#'   \code{formula}.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}.
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using \code{document_weights}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{TS_fit}-class list containing the following elements, many of
#'   which are hidden for \code{print}ing, but are accessible:
#'   \describe{
#'     \item{data}{\code{data} input to the function.}
#'     \item{formula}{\code{\link[stats]{formula}} input to the function.}
#'     \item{nchangepoints}{\code{nchangepoints} input to the function.}
#'     \item{weights}{\code{weights} input to the function.}
#'     \item{control}{\code{control} input to the function.}
#'     \item{lls}{Iteration-by-iteration 
#'                \link[=logLik.multinom_TS_fit]{logLik} values for the
#'                 full time series fit by \code{\link{multinom_TS}}.}
#'     \item{rhos}{Iteration-by-iteration change point estimates from
#'                 \code{\link{est_changepoints}}.}
#'     \item{etas}{Iteration-by-iteration marginal regressor estimates from
#'                 \code{\link{est_regressors}}, which have been 
#'                 unconditioned with respect to the change point locations.}
#'     \item{ptMCMC_diagnostics}{ptMCMC diagnostics, 
#'                                see \code{\link{diagnose_ptMCMC}}}
#'     \item{rho_summary}{Summary table describing \code{rhos} (the change
#'                        point locations), 
#'                        see \code{\link{summarize_rhos}}.}
#'     \item{rho_vcov}{Variance-covariance matrix for the estimates of
#'                      \code{rhos} (the change point locations), see 
#'                      \code{\link{measure_rho_vcov}}.}
#'     \item{eta_summary}{Summary table describing \code{ets} (the 
#'                        regressors), 
#'                        see \code{\link{summarize_etas}}.}
#'     \item{eta_vcov}{Variance-covariance matrix for the estimates of
#'                      \code{etas} (the regressors), see 
#'                      \code{\link{measure_eta_vcov}}.}
#'     \item{logLik}{Across-iteration average of log-likelihoods 
#'                    (\code{lls}).}
#'     \item{nparams}{Total number of parameters in the full model,
#'                    including the change point locations and regressors.}
#'     \item{deviance}{Penalized negative log-likelihood, based on 
#'                     \code{logLik} and \code{nparams}.}
#'   }
#'
#' @references
#'   Bishop, C. M. 2006. \emph{Pattern Recognition and Machine Learning}. 
#'    Springer, New York, NY, USA.
#'
#'   Christensen, E., D. J. Harris, and S. K. M. Ernest. 2018.
#'   Long-term community change through multiple rapid transitions in a 
#'   desert rodent community. \emph{Ecology} \strong{99}:1523-1529. 
#'   \href{https://doi.org/10.1002/ecy.2373}{link}.
#'
#'   Earl, D. J. and M. W. Deem. 2005. Parallel tempering: theory, 
#'   applications, and new perspectives. \emph{Physical Chemistry Chemical 
#'   Physics} \strong{7}: 3910-3916.
#'   \href{https://doi.org/10.1039/B509983H}{link}.
#'
#'   McCullagh, P. and J. A. Nelder. 1989. \emph{Generalized Linear Models}.
#'   2nd Edition. Chapman and Hall, New York, NY, USA.
#'
#'   Ripley, B. D. 1996. \emph{Pattern Recognition and Neural Networks}. 
#'   Cambridge University Press, Cambridge, UK.
#'
#'   Ruggieri, E. 2013. A Bayesian approach to detecting change points in 
#'   climactic records. \emph{International Journal of Climatology}
#'   \strong{33}:520-528.
#'   \href{https://doi.org/10.1002/joc.3447}{link}.
#'
#'   Venables, W. N. and B. D. Ripley. 2002. \emph{Modern and Applied
#'   Statistics with S}. Fourth Edition. Springer, New York, NY, USA.
#'
#'   Western, B. and M. Kleykamp. 2004. A Bayesian change point model for 
#'   historical time series analysis. \emph{Political Analysis}
#'   \strong{12}:354-374.
#'   \href{https://doi.org/10.1093/pan/mph023}{link}.
#'
#' @examples 
#' \dontrun{
#'   data(rodents)
#'   document_term_table <- rodents$document_term_table
#'   document_covariate_table <- rodents$document_covariate_table
#'   LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
#'   data <- document_covariate_table
#'   data$gamma <- LDA_models@gamma
#'   weights <- document_weights(document_term_table)
#'   TSmod <- TS(data, gamma ~ 1, nchangepoints = 1, "newmoon", weights)
#' }
#'
#' @export
#'
TS <- function(data, formula, nchangepoints, timename = "time",  
               weights = NULL, control = list()){
  check_TS_inputs(data, formula, nchangepoints, timename, weights, control)
  control <- do.call("TS_control", control)
  set.seed(control$seed)
  data <- data[order(data[,timename]), ]
  rho_dist <- est_changepoints(data, formula, nchangepoints, timename, 
                               weights, control)
  eta_dist <- est_regressors(rho_dist, data, formula, timename, weights, 
                             control)
  package_TS(data, formula, timename, weights, control, rho_dist, eta_dist)
}

#' @rdname TS
#'
#' @description \code{check_TS_inputs} checks that the inputs to 
#'   \code{TS} are of proper classes for a full analysis.
#' 
#' @export
#'
check_TS_inputs <- function(data, formula, nchangepoints, timename, weights, 
                            control){
  check_formula(data, formula)  
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_timename(data, timename)
  check_control(control)
}

#' @title Summarize the Time Series model 
#'
#' @description Calculate relevant summaries for the run of a Time Series
#'   model within \code{\link{TS}} and package the output as a
#'   \code{TS_fit}-class object.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. 
#'
#' @param formula \code{\link[stats]{formula}} defining the regression between
#'   relationship the change points. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (multinomial) response variable must be a set of
#'   columns in \code{data}, as verified by \code{\link{check_formula}}.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. 
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using \code{document_weights}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @param rho_dist List of saved data objects from the ptMCMC estimation of
#'   change point locations returned by \code{\link{est_changepoints}}
#'   (unless \code{nchangepoints} is 0, then \code{NULL}).
#'
#' @param eta_dist Matrix of draws (rows) from the marginal posteriors of the 
#'   coefficients across the segments (columns), as estimated by
#'   \code{\link{est_regressors}}. 
#'
#' @return \code{TS_fit}-class list containing the following elements, many of
#'   which are hidden for \code{print}ing, but are accessible:
#'   \describe{
#'     \item{data}{\code{data} input to the function.}
#'     \item{formula}{\code{\link[stats]{formula}} input to the function.}
#'     \item{nchangepoints}{\code{nchangepoints} input to the function.}
#'     \item{weights}{\code{weights} input to the function.}
#'     \item{timename}{\code{timename} input to the function.}
#'     \item{control}{\code{control} input to the function.}
#'     \item{lls}{Iteration-by-iteration 
#'                \link[=logLik.multinom_TS_fit]{logLik} values for the
#'                 full time series fit by \code{\link{multinom_TS}}.}
#'     \item{rhos}{Iteration-by-iteration change point estimates from
#'                 \code{\link{est_changepoints}}.}
#'     \item{etas}{Iteration-by-iteration marginal regressor estimates from
#'                 \code{\link{est_regressors}}, which have been 
#'                 unconditioned with respect to the change point locations.}
#'     \item{ptMCMC_diagnostics}{ptMCMC diagnostics, 
#'                                see \code{\link{diagnose_ptMCMC}}}
#'     \item{rho_summary}{Summary table describing \code{rhos} (the change
#'                        point locations), 
#'                        see \code{\link{summarize_rhos}}.}
#'     \item{rho_vcov}{Variance-covariance matrix for the estimates of
#'                      \code{rhos} (the change point locations), see 
#'                      \code{\link{measure_rho_vcov}}.}
#'     \item{eta_summary}{Summary table describing \code{ets} (the 
#'                        regressors), 
#'                        see \code{\link{summarize_etas}}.}
#'     \item{eta_vcov}{Variance-covariance matrix for the estimates of
#'                      \code{etas} (the regressors), see 
#'                      \code{\link{measure_eta_vcov}}.}
#'     \item{logLik}{Across-iteration average of log-likelihoods 
#'                    (\code{lls}).}
#'     \item{nparams}{Total number of parameters in the full model,
#'                    including the change point locations and regressors.}
#'     \item{deviance}{Penalized negative log-likelihood, based on 
#'                     \code{logLik} and \code{nparams}.}
#'   }
#'
#' @export
#'
package_TS <- function(data, formula, timename, weights, control, rho_dist, 
                         eta_dist){

  check_formula(data, formula) 
  check_weights(weights)
  check_control(control)
  check_timename(data, timename)
  control <- do.call("TS_control", control)
  nchangepoints <- dim(rho_dist$cpts)[1]
  if (is.null(nchangepoints)){
    nchangepoints <- 0
    mod <- multinom_TS(data, formula, changepoints = NULL, timename, weights,
                       control)
    mod <- mod[[1]][[1]]
    lls <- as.numeric(logLik(mod))
    rhos <- NULL
  } else{
    lls <- rho_dist$lls[1, ]
    rhos <- t(array(rho_dist$cpts[ , 1, ], dim = dim(rho_dist$cpts)[c(1, 3)]))
  }

  ptMCMC_diagnostics <- diagnose_ptMCMC(rho_dist)
  rho_summary <- summarize_rhos(rhos, control)
  rho_vcov <- measure_rho_vcov(rhos)
  eta_summary <- summarize_etas(eta_dist, control)
  eta_vcov <- measure_eta_vcov(eta_dist)

  logLik <- mean(lls)
  ncoefs <- ncol(eta_dist)
  nparams <- nchangepoints + ncoefs 
  deviance <- -2 * logLik + 2 * nparams

  out <- list(data = data, formula = formula, nchangepoints = nchangepoints,
              timename = timename, weights = weights,
              control = control, lls = lls, rhos = rhos,
              etas = eta_dist, ptMCMC_diagnostics = ptMCMC_diagnostics,
              rho_summary = rho_summary, rho_vcov = rho_vcov,
              eta_summary = eta_summary, eta_vcov = eta_vcov,
              logLik = logLik, nparams = nparams, deviance = deviance)
  class(out) <- c("TS_fit", "list")
  to_hide <- c("data", "weights", "control", "lls", "rhos", "etas", 
               "rho_vcov", "eta_vcov")
  if (nchangepoints == 0){
    to_hide <- c(to_hide, "ptMCMC_diagnostics", "rho_summary")
  }
  attr(out, "hidden") <- to_hide
  out
}

#' @title Print a Time Series model fit
#'
#' @description Convenience function to print only the most important 
#'   components of a \code{TS_fit}-class object fit by 
#'   \code{\link{TS}}.
#'
#' @param x Class \code{TS_fit} object to be printed.
#'
#' @param ... Not used, simply included to maintain method compatability.
#'
#' @export
#'
print.TS_fit <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !names(x) %in% hid
  print(x[notHid])
}

#' @title Summarize the regressor (eta) distributions
#'
#' @description \code{summarize_etas} calculates summary statistics for each
#'   of the chunk-level regressors. 
#'   \cr \cr
#'   \code{measure_ets_vcov} generates the variance-covariance matrix for 
#'   the regressors.
#'
#' @param etas Matrix of regressors (columns) across iterations of the 
#'   ptMCMC (rows), as returned from \code{\link{est_regressors}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{summarize_etas}: table of summary statistics for chunk-level
#'   regressors including mean, median, mode, posterior interval, standard
#'   deviation, MCMC error, autocorrelation, and effective sample size for 
#'   each regressor.
#' 
#' @export 
#'
summarize_etas <- function(etas, control = list()){
  check_control(control)
  control <- do.call("TS_control", control)
  if (!is.matrix(etas)){
    stop("etas should be a matrix")
  }
  prob <- control$summary_prob
  Mean <- round(apply(etas, 2, mean), 4)
  Median <- round(apply(etas, 2, median), 4)
  SD <- round(apply(etas, 2, sd), 4)
  MCMCerr <- round(SD / sqrt(nrow(etas)), 4)
  HPD <- HPDinterval(as.mcmc(etas), prob = prob)
  Lower <- round(HPD[ , "lower"], 4)
  Upper <- round(HPD[ , "upper"], 4)
  AC10 <- tryCatch(t(round(autocorr.diag(as.mcmc(etas), lag = 10), 4)),
                   error = function(x){"-"})
  ESS <- effectiveSize(etas)
  out <- data.frame(Mean, Median, Lower, Upper, SD, MCMCerr, AC10, ESS)
  colnames(out)[3:4] <- paste0(c("Lower_", "Upper_"), paste0(prob*100, "%"))
  colnames(out)[7] <- "AC10"
  rownames(out) <- colnames(etas)
  out
}

#' @rdname summarize_etas 
#'
#' @return \code{measure_eta_vcov}: variance-covariance matrix for chunk-level
#'   regressors.
#'
#' @export
#'
measure_eta_vcov <- function(etas){
  if (!is.matrix(etas)){
    stop("expecting etas to be a matrix")
  }
  out <- var(etas)
  colnames(out) <- colnames(etas)
  rownames(out) <- colnames(etas)
  out
}

#' @title Summarize the rho distributions
#'
#' @description \code{summarize_rho} calculates summary statistics for each
#'   of the change point locations.
#'   \cr \cr
#'   \code{measure_rho_vcov} generates the variance-covariance matrix for the 
#'   change point locations.
#'
#' @param rhos Matrix of change point locations (columns) across iterations of 
#'   the ptMCMC (rows) or \code{NULL} if no change points are in the model,
#'   as returned from \code{\link{est_changepoints}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{summarize_rhos}: table of summary statistics for change point
#'   locations including mean, median, mode, posterior interval, standard
#'   deviation, MCMC error, autocorrelation, and effective sample size for 
#'   each change point location.
#' 
#' @export 
#'
summarize_rhos <- function(rhos, control = list()){
  check_control(control)
  control <- do.call("TS_control", control)
  if (is.null(rhos)) {
    return()
  }
  prob <- control$summary_prob
  Mean <- round(apply(rhos, 2, mean), 2)
  Median <- apply(rhos, 2, median)
  Mode <- apply(rhos, 2, modalvalue)
  SD <- round(apply(rhos, 2, sd), 2)
  MCMCerr <- round(SD / sqrt(nrow(rhos)), 4)
  HPD <- HPDinterval(as.mcmc(rhos), prob = prob)
  Lower <- HPD[ , "lower"]
  Upper <- HPD[ , "upper"]
  AC10 <- t(round(autocorr.diag(as.mcmc(rhos), lag = 10), 4))
  ESS <- effectiveSize(rhos)
  out <- data.frame(Mean, Median, Mode, Lower, Upper, SD, MCMCerr, AC10, ESS)
  colnames(out)[4:5] <- paste0(c("Lower_", "Upper_"), paste0(prob*100, "%"))
  colnames(out)[8] <- "AC10"
  rownames(out) <- sprintf("Changepoint_%d", seq_len(nrow(out)))
  out
}

#' @rdname summarize_rhos 
#'
#' @return \code{measure_rho_vcov}: variance-covariance matrix for change 
#'   point locations.
#'
#' @export
#'
measure_rho_vcov <- function(rhos){
  if (is.null(rhos)) {
    return()
  }
  if (!is.matrix(rhos)){
    stop("expecting rhos to be a matrix")
  }
  out <- var(rhos)
  colnames(out) <- sprintf("CP_%d", 1:dim(out)[1])
  rownames(out) <- sprintf("CP_%d", 1:dim(out)[2])
  out
}

#' @title Estimate the distribution of regressors, unconditional on the
#'   change point locations
#'
#' @description This function uses the marginal posterior distributions of
#'   the change point locations (estimated by \code{\link{est_changepoints}})
#'   in combination with the conditional (on the change point locations) 
#'   posterior distributions of the regressors (estimated by
#'   \code{\link{multinom_TS}}) to estimate the marginal posterior 
#'   distribution of the regressors, unconditional on the change point 
#'   locations.
#'
#' @details The general approach follows that of Western and Kleykamp
#'   (2004), although we note some important differences. Our regression
#'   models are fit indpendently for each chunk (segment of time), and 
#'   therefore the variance-covariance matrix for the full model 
#'   has \code{0} entries for covariances between regressors in different
#'   chunks of the time series. Further, because the regression model here
#'   is a standard (non-hierarchical) softmax (Ripley 1996, Venables and 
#'   Ripley 2002, Bishop 2006), there is no error term in the regression  
#'   (as there is in the normal model used by Western and Kleykamp 2004), 
#'   and so the posterior distribution used here is a multivariate normal,
#'   as opposed to a multivariate t, as used by Western and Kleykamp (2004).
#'
#' @param rho_dist List of saved data objects from the ptMCMC estimation of
#'   change point locations (unless \code{nchangepoints} is 0, then 
#'   \code{NULL}) returned from \code{\link{est_changepoints}}.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. 
#'
#' @param formula \code{\link[stats]{formula}} defining the regression between
#'   relationship the change points. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (multinomial) response variable must be a set of
#'   columns in \code{data}, as verified by \code{\link{check_formula}}.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. 
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using \code{document_weights}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return Matrix of draws (rows) from the marginal posteriors of the 
#'   coefficients across the segments (columns). 
#'
#' @references
#'   Bishop, C. M. 2006. \emph{Pattern Recognition and Machine Learning}. 
#'    Springer, New York, NY, USA.
#'
#'   Ripley, B. D. 1996. \emph{Pattern Recognition and Neural Networks}. 
#'   Cambridge University Press, Cambridge, UK.
#'
#'   Venables, W. N. and B. D. Ripley. 2002. \emph{Modern and Applied
#'   Statistics with S}. Fourth Edition. Springer, New York, NY, USA.
#'
#'   Western, B. and M. Kleykamp. 2004. A Bayesian change point model for 
#'   historical time series analysis. \emph{Political Analysis}
#'   \strong{12}:354-374.
#'   \href{https://doi.org/10.1093/pan/mph023}{link}.
#'
#' @export
#'
est_regressors <- function(rho_dist, data, formula, timename, weights, 
                           control = list()){
  check_formula(data, formula)  
  check_weights(weights)
  check_control(control)
  control <- do.call("TS_control", control)
  if (!is.null(rho_dist)){
    if (any(names(rho_dist)[1:3] != c("cpts", "lls", "ids"))){
      stop("expecting rho_dist to have elements cpts, lls, ids")
    }
  }
  if (is.null(rho_dist)){
    mod <- multinom_TS(data, formula, changepoints = NULL, timename, weights, 
                       control)
    mod <- mod[[1]][[1]]
    mv <- as.vector(t(coef(mod)))
    vcv <- mirror_vcov(mod)
    eta <- rmvnorm(control$nit, mv, vcv)
    seg_names <- rep(1, ncol(vcv))
    coef_names <- colnames(vcv)
    colnames(eta) <- paste(seg_names, coef_names, sep = "_")
    return(eta)
  }

  focal_rho <- rho_dist$cpts[ , 1, ]
  nchangepts <- dim(rho_dist$cpts)[1]
  if (nchangepts == 1){
    collapsedrho <- focal_rho
  } else{
    collapsedrho <- apply(focal_rho, 2, paste, collapse = "_")
  }
  freq_r <- table(collapsedrho)
  unique_r <- names(freq_r)
  nr <- length(unique_r)
  n_topic <- ncol(data$gamma)
  n_covar <- length(attr(terms(formula), "term.labels"))
  n_eta_segment <- (n_topic - 1) * (n_covar + 1)
  n_changept <- dim(rho_dist$cpts)[1]
  n_segment <- n_changept + 1
  n_eta <- n_eta_segment * n_segment 
  eta <- matrix(NA, nrow = control$nit, ncol = n_eta)
  pbar <- prep_pbar(control, "eta", nr)

  for(i in 1:nr){
    update_pbar(pbar, control)
    cpts <- as.numeric(strsplit(unique_r[i], "_")[[1]])
    mods <- multinom_TS(data, formula, cpts, timename, weights, control)
    ndraws <- freq_r[i]
    colindex1 <- 1
    for(j in 1:n_segment){
      colindex2 <- colindex1 + n_eta_segment - 1
      seg_mod <- mods[[1]][[j]]
      mv <- as.vector(t(coef(seg_mod)))
      vcv <- mirror_vcov(seg_mod)
      drawn <- rmvnorm(ndraws, mv, vcv)    
      rows_in <- which(collapsedrho == unique_r[i])
      cols_in <- colindex1:colindex2
      eta[rows_in, cols_in] <- drawn
      colindex1 <- colindex2 + 1
    }
  }
  seg_names <- rep(1:n_segment, each = n_eta_segment)
  coef_names <- rep(colnames(vcv), n_segment)
  colnames(eta) <- paste(seg_names, coef_names, sep = "_")
  eta
}


#' @title Use ptMCMC to estimate the distribution of change point locations
#'
#' @description This function executes ptMCMC-based estimation of the 
#'   change point location distributions for multinomial Time Series analyses.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. 
#'
#' @param formula \code{\link[stats]{formula}} defining the regression between
#'   relationship the change points. Any 
#'   predictor variable included must also be a column in 
#'   \code{data} and any (multinomial) response variable must be a set of
#'   columns in \code{data}, as verified by \code{\link{check_formula}}.
#'
#' @param nchangepoints \code{integer} corresponding to the number of 
#'   change points to include in the model. 0 is a valid input (corresponding
#'   to no change points, so a singular time series model), and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segementation of the 
#'   time series into chunks fit with separate models dictated by 
#'   \code{formula}.
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. 
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{multinom_TS} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using \code{document_weights}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return List of saved data objects from the ptMCMC estimation of
#'   change point locations (unless \code{nchangepoints} is 0, then 
#'   \code{NULL} is returned).
#'
#' @export
#'
est_changepoints <- function(data, formula, nchangepoints, timename, weights, 
                             control = list()){
  check_TS_inputs(data, formula, nchangepoints, timename, weights, control)
  control <- do.call("TS_control", control)
  if (nchangepoints == 0){
    return(NULL)
  }
  saves <- prep_saves(nchangepoints, control)
  inputs <- prep_ptMCMC_inputs(data, formula, nchangepoints, timename, 
                               weights, control)
  cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
  ids <- prep_ids(control)
  pbar <- prep_pbar(control, "rho")

  for(i in 1:control$nit){
    update_pbar(pbar, control)
    steps <- step_chains(i, cpts, inputs)
    swaps <- swap_chains(steps, inputs, ids)
    saves <- update_saves(i, saves, steps, swaps)
    cpts <- update_cpts(cpts, swaps)
    ids <- update_ids(ids, swaps)
  }
  process_saves(saves, control)
}

#' @title Initialize and tick through the progress bar
#'
#' @description \code{prep_pbar} creates and \code{update_pbar} steps
#'   through the progress bars (if desired) in \code{\link{TS}}
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model including the parallel tempering Markov Chain 
#'   Monte Carlo (ptMCMC) controls. Values not input assume defaults set by 
#'   \code{\link{TS_control}}. Of use here is \code{quiet} which is a
#'   a \code{logical} indicator of whether there should be information 
#'   (i.e. the progress bar) printed during the run or not. Default is 
#'   \code{TRUE}.
#'
#' @param bar_type "rho" (for change point locations) or "eta" (for 
#'   regressors).
#'
#' @param nr \code{integer} number of unique realizations, needed when
#'   \code{bar_type} = "eta".
#'
#' @return \code{prep_pbar}: the initialized progress bar object.
#'
#' @export
#'
prep_pbar <- function(control = list(), bar_type = "rho", 
                      nr = NULL){
  check_control(control)
  control <- do.call("TS_control", control)
  if (control$quiet){
    return()
  }
  if (!(bar_type %in% c("eta", "rho"))){
    stop("bar_type must be eta or rho")
  }
  if (!is.null(nr)){
    if (!is.numeric(nr) || any(nr %% 1 != 0)){
      stop("nr must be integer-valued")
    }
  }
  form <- "  [:bar] :percent eta: :eta"
  if (bar_type == "rho"){
    cat("  Estimating changepoint distribution \n") 
    out <- progress_bar$new(form, control$nit, width = 60)
  }
  if (bar_type == "eta"){
    cat("  Estimating regressor distribution \n") 
    out <- progress_bar$new(form, nr, width = 60)
  }
  out
}

#' @rdname prep_pbar
#'
#' @param pbar The progress bar object returned from \code{prep_pbar}.
#'
#' @return \code{update_pbar}: the ticked-forward \code{pbar}.
#'
#' @export
#'
update_pbar <- function(pbar, control = list()){
  if (!("progress_bar" %in% class(pbar))){
    stop("pbar must be of class progress_bar")
  }
  check_control(control)
  control <- do.call("TS_control", control)
  if (control$quiet){
    return()
  }
  pbar$tick()
}


#' @title Check that a formula is proper
#' 
#' @description Check that \code{formula} is actually a 
#'   \code{\link[stats]{formula}} and that the
#'   response and predictor variabless are all included in \code{data}.
#'   
#' @param formula \code{formula} to evaluate.
#'
#' @param data \code{data.frame} including [1] the time variable (indicated 
#'   in \code{timename}), [2] the predictor variables (required by
#'   \code{formula}) and [3], the multinomial response variable (indicated in
#'   \code{formula}) as verified by \code{\link{check_timename}} and 
#'   \code{\link{check_formula}}. Note that the response variables should be
#'   formatted as a \code{data.frame} object named as indicated by the 
#'   \code{response} entry in the \code{control} list, such as \code{gamma} 
#'   for a standard TS analysis on LDA output. 
#' 
#' @export
#'
check_formula <- function(data, formula){

  if (!is(formula, "formula")){
    stop("formula does not contain a single formula")
  }

  respLoc <- attr(terms(formula), "response")
  if (respLoc == 0){
    stop("formula inputs should include response variable")
  }  

  resp <- as.character(attr(terms(formula), "variables"))[-1][respLoc]
  pred <- attr(terms(formula), "term.labels")
  if (!resp %in% colnames(data)){
    stop("formula includes response not present in data")
  }
  if (!all(pred %in% colnames(data))){
    misses <- pred[which(pred %in% colnames(data) == FALSE)]
    mis <- paste(misses, collapse = ", ")
    stop(paste0("formula includes predictors not present in data: ", mis))
  }
}

#' @title Create the controls list for the Time Series model
#'
#' @description This function provides a simple creation and definition of a
#'   list used to control the time series model fit occurring within 
#'   \code{\link{TS}}. 
#'
#' @param memoise \code{logical} indicator of whether the multinomial 
#'   functions should be memoised (via \code{\link[memoise]{memoise}}). 
#'   Memoisation happens to both \code{\link{multinom_TS}} and 
#'   \code{\link{multinom_TS_chunk}}.
#'
#' @param response \code{character} element indicating the response variable 
#'   used in the time series. 
#'
#' @param lambda \code{numeric} "weight" decay term used to set the prior
#'   on the regressors within each chunk-level model. Defaults to 0, 
#'   corresponding to a fully vague prior.
#'
#' @param measurer,selector Function names for use in evaluation of the TS
#'   models. \code{measurer} is used to create a value for each model
#'   and \code{selector} operates on the values to choose the model. 
#'
#' @param ntemps \code{integer} number of temperatures (chains) to use in the 
#'   ptMCMC algorithm.
#'
#' @param penultimate_temp Penultimate temperature in the ptMCMC sequence.
#'
#' @param ultimate_temp Ultimate temperature in the ptMCMC sequence.
#'
#' @param q Exponent controlling the ptMCMC temperature sequence from the 
#'   focal chain (reference with temperature = 1) to the penultimate chain. 0
#'   (default) implies a geometric sequence. 1 implies squaring before 
#'   exponentiating.
#'
#' @param nit \code{integer} number of iterations (steps) used in the ptMCMC
#'   algorithm.
#'
#' @param magnitude Average magnitude (defining a geometric distribution)
#'   for the proposed step size in the ptMCMC algorithm.
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param burnin \code{integer} number of iterations to remove from the 
#'   beginning of the ptMCMC algorithm.
#'
#' @param thin_frac Fraction of iterations to retain, must be \eqn{(0, 1]},
#'   and the default value of 1 represents no thinning.
#'
#' @param summary_prob Probability used for summarizing the posterior 
#'   distributions (via the highest posterior density interval, see
#'   \code{\link[coda]{HPDinterval}}).
#'
#' @param seed Input to \code{set.seed} for replication purposes.
#'
#' @return \code{list}, with named elements corresponding to the arguments.
#'
#' @export
#'
TS_control <- function(memoise = TRUE, response = "gamma", lambda = 0, 
                       measurer = AIC, selector = min, ntemps = 6, 
                       penultimate_temp = 2^6, ultimate_temp = 1e10, q = 0, 
                       nit = 1e4, magnitude = 12, quiet = FALSE, burnin = 0, 
                       thin_frac = 1, summary_prob = 0.95, seed = NULL){
  list(memoise = memoise, response = response, lambda = lambda, 
       measurer = measurer, selector = selector, ntemps = ntemps, 
       penultimate_temp = penultimate_temp, ultimate_temp = ultimate_temp, 
       q = q, nit = nit, magnitude = magnitude, quiet = quiet, 
       burnin = burnin, thin_frac = thin_frac, summary_prob = summary_prob,
       seed = seed)

}

#' @title Determine the AIC (deviance) value of a Time Series model
#'
#' @description Convenience function to extract the AIC (deviance) element of 
#'   \code{TS_fit}-class object fit by \code{\link{multinom_TS}}.
#'
#' @param object Class \code{TS_fit} object to be evaluated.
#'
#' @param ... Not used, simply included to maintain method compatability.
#'
#' @param k Not used, simply included to maintain method compatability.
#'
#' @export
#'
AIC.TS_fit <- function(object, ..., k = 2){
  object$deviance
}
