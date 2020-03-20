#' @title Estimate a Time Series model sequentially
#'
#' @description This set of functions estimates the Time Series model
#'   by sequential methods that first estimate the change point locations
#'   with full flexibility of the regressor models between change points,
#'   then estimate the regressors between the change points, intially 
#'   conditional on their locations, but with marginal estimation to produce
#'   regressor values unconditional on change point locations. 
#'   \code{sequential_TS} combines each stage of the model estimation and 
#'     packages the model results in a consistent output. 
#'   \code{est_changepoints} estimates the change point location 
#'     distributions for multinomial Time Series analyses.
#'   \code{est_regressors} uses the marginal posterior distributions of
#'     the change point locations (estimated by 
#'     \code{\link{est_changepoints}}) in combination with the conditional 
#'     (on the change point locations) posterior distributions of the 
#'     regressors (estimated by a \code{<response>_TS} function) to 
#'     estimate the marginal posterior distribution of the regressors, 
#'     unconditional on the change point locations.
#'   \code{package_sequential_TS} calculates relevant summaries for the run of
#'     a sequenial Time Series model within \code{\link{sequential_TS}} and 
#'     packages the output as a \code{TS}-class object.  
#'   \code{sequential_TS_msg} produces a specific message about the model
#'     being run. 
#'
#' @param rho_dist \code{list} of saved data objects from the estimation of
#'   change point locations (unless \code{nchangepoints} is 0, then 
#'   \code{NULL}) returned from \code{\link{est_changepoints}}.
#'
#' @param eta_dist \code{matrix} of draws (rows) from the marginal posteriors
#'   of the coefficients across the segments (columns), as estimated by
#'   \code{\link{est_regressors}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @details The general approach follows that of Western and Kleykamp
#'   (2004), although we note some important differences. Our regression
#'   models are fit independently for each chunk (segment of time), and 
#'   therefore the variance-covariance matrix for the full model 
#'   has \code{0} entries for covariances between regressors in different
#'   chunks of the time series. Further, because the regression model here
#'   is a standard (non-hierarchical) softmax (Ripley 1996, Venables and 
#'   Ripley 2002, Bishop 2006), there is no error term in the regression  
#'   (as there is in the normal model used by Western and Kleykamp 2004), 
#'   and so the posterior distribution used here is a multivariate normal,
#'   as opposed to a multivariate t, as used by Western and Kleykamp (2004).
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
#' @return 
#'   \code{est_changepoints}: \code{list} of saved data objects from the 
#'     estimation of change point locations, uunless \code{nchangepoints} 
#'     is 0, then \code{NULL}.
#'   \code{est_regressors}: \code{matrix} of draws (rows) from the marginal 
#'     posteriors of the coefficients across the segments (columns). 
#'   \code{sequential_TS} and \code{package_sequential_TS}: 
#'     \code{TS}-class list containing the following elements, many of
#'      which are hidden for \code{print}ing, but are accessible:
#'     \describe{
#'       \item{data}{\code{data} input to the function.}
#'       \item{formula}{\code{\link[stats]{formula}} input to the function.}
#'       \item{nchangepoints}{\code{nchangepoints} input to the function.}
#'       \item{weights}{\code{weights} input to the function.}
#'       \item{timename}{\code{timename} input to the function.}
#'       \item{control}{\code{control} input to the function.}
#'       \item{lls}{Iteration-by-iteration 
#'                  \link[=logLik.multinom_TS_fit]{logLik} values for the
#'                   full time series fit by \code{\link{multinom_TS}}.}
#'       \item{rhos}{Iteration-by-iteration change point estimates from
#'                   \code{\link{est_changepoints}} and diagnostics.}
#'       \item{focal_rhos}{Simplified object of just the change point 
#'                         locations of interest.}
#'       \item{etas}{Iteration-by-iteration marginal regressor estimates from
#'                   \code{\link{est_regressors}}, which have been 
#'                   unconditioned with respect to change point locations.}
#'       \item{rho_summary}{Summary table describing \code{rhos} (the change
#'                          point locations), see 
#'                          \code{\link{summarize_rhos}}.}
#'       \item{rho_vcov}{Variance-covariance matrix for the estimates of
#'                      \code{rhos} (the change point locations), see 
#'                      \code{\link{measure_rho_vcov}}.}
#'       \item{eta_summary}{Summary table describing \code{ets} (the 
#'                          regressors), see 
#'                          \code{\link{summarize_etas}}.}
#'       \item{eta_vcov}{Variance-covariance matrix for the estimates of
#'                      \code{etas} (the regressors), see 
#'                      \code{\link{measure_eta_vcov}}.}
#'       \item{logLik}{Across-iteration average of log-likelihoods 
#'                    (\code{lls}).}
#'       \item{nparams}{Total number of parameters in the full model,
#'                      including the change point locations and regressors.}
#'     }
#'
#' @export
#'
sequential_TS <- function(TS, control = list()){
  control <- do.call("TS_control", control)
  sequential_TS_msg(TS = TS, control = control)
  rho_dist <- est_changepoints(TS = TS, control = control)
  eta_dist <- est_regressors(rho_dist = rho_dist, TS = TS, control = control)
  package_sequential_TS(TS = TS, rho_dist = rho_dist, eta_dist = eta_dist, 
                        control = control)
}

#'
#' @rdname sequential_TS
#'
#' @export
#'
package_sequential_TS <- function(TS, rho_dist, eta_dist, control = list()){


  if(is.null(rho_dist)){
    focal_rho_dist <- NULL
    data <- TS$data$train$ts_data
    fun <- eval(parse(text = paste0(TS$response, "_TS")))
    args <- list(data = data, formula = TS$formula, changepoints = NULL, 
                 timename = TS$timename, weights = TS$weights, 
                 control = control$method_args)
    mod <- soft_call(fun, args, TRUE)
    lls <- as.numeric(logLik(mod))



  } else{
    vals <- rho_dist$cpts[ , 1, , drop = FALSE]
    dims <- dim(rho_dist$cpts)[c(1, 3)]
    lls <- rho_dist$lls[1, ]
    focal_rho_dist <- t(array(vals, dim = dims))
  }


  rho_summary <- summarize_rhos(rhos = focal_rho_dist, control = control)
  rho_vcov <- measure_rho_vcov(rhos = focal_rho_dist)
  eta_summary <- summarize_etas(etas = eta_dist, control = control)
  eta_vcov <- measure_eta_vcov(etas = eta_dist)

  logLik <- mean(lls)
  ncoefs <- ncol(eta_dist)
  nparams <- TS$nchangepoints + ncoefs 

  out <- update_list(TS, focal_rhos = focal_rho_dist, rhos = rho_dist,
                     etas = eta_dist, rho_summary = rho_summary,
                     rho_vcov = rho_vcov, eta_summary = eta_summary,
                     eta_vcov = eta_vcov, logLik = logLik, nparams = nparams)
  class(out) <- c("TS", "list")
  to_hide <- c("data", "weights", "control", "lls", "rhos", "etas", 
               "focal_rhos", "rho_vcov", "eta_vcov")
  if (TS$nchangepoints == 0){
    to_hide <- c(to_hide, "rho_summary")
  }
  attr(out, "hidden") <- to_hide
  out  
}

#'
#' @rdname sequential_TS
#'
#' @export
#'
est_changepoints <- function(TS, control = list()){
  if (TS$nchangepoints == 0){
    return(NULL)
  }
  fun <- control$method
  args <- list(TS = TS, control = control$method_args)
  soft_call(fun = fun, args = args, soften = control$soften)
}

#'
#' @rdname sequential_TS
#'
#' @export
#'
est_regressors <- function(rho_dist, TS, control = list()){
  data <- TS$data$train$ts_data
  if(is.null(rho_dist)){
    fun <- eval(parse(text = paste0(TS$response, "_TS")))
    args <- list(data = data, formula = TS$formula, changepoints = NULL, 
                 timename = TS$timename, weights = TS$weights, 
                 control = control$method_args)
    mod <- soft_call(fun, args, TRUE)

    mod <- mod[[1]][[1]]
    mv <- as.vector(t(coef(mod)))
    vcv <- mirror_vcov(mod)
    eta <- rmvnorm(control$method_args$nit, mv, vcv)
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
  n_covar <- length(attr(terms(TS$formula), "term.labels"))
  n_eta_segment <- (n_topic - 1) * (n_covar + 1)
  n_changept <- dim(rho_dist$cpts)[1]
  n_segment <- n_changept + 1
  n_eta <- n_eta_segment * n_segment 
  n_iter <- dim(rho_dist$cpts)[3]
  eta <- matrix(NA, nrow = n_iter, ncol = n_eta)
  pbar <- prep_pbar(control = control, bar_type = "eta", nr = nr)

  for(i in 1:nr){
    update_pbar(pbar = pbar, control = control)
    cpts <- as.numeric(strsplit(unique_r[i], "_")[[1]])


    data <- TS$data$train$ts_data
    fun <- eval(parse(text = paste0(TS$response, "_TS")))
    fun <- memoise_fun(fun, control$method_args$memoise)
    args <- list(data = data, formula = TS$formula, changepoints = cpts, 
                 timename = TS$timename, weights = TS$weights, 
                 control = control$method_args)
    mod <- soft_call(fun, args, TRUE)


    ndraws <- freq_r[i]
    colindex1 <- 1
    for(j in 1:n_segment){
      colindex2 <- colindex1 + n_eta_segment - 1
      seg_mod <- mod[[1]][[j]]
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

#'
#' @rdname sequential_TS
#'
#' @export
#'
sequential_TS_msg <- function(TS, control = list()){
  subset_msg <- paste0("  - data subset ", TS$data_subset)
  topic_msg <- paste0(", ", TS$topics, " topics")
  rep_msg <- paste0(", replicate ", TS$rep)

  formula_msg <- paste0(", ", deparse(TS$formula))
  nchangepoints <- TS$nchangepoints
  txt <- ifelse(nchangepoints == 1, " change point", " change points")
  changepoints_msg <- paste0(", ", nchangepoints, txt)
  msg <- paste0(subset_msg, topic_msg, rep_msg, formula_msg, changepoints_msg)
  messageq(msg, control$quiet)
}


#' @title Summarize the regressor (eta) distributions of a time series model
#'
#' @description \code{summarize_etas} calculates summary statistics for each
#'   of the chunk-level regressors. 
#'   \cr \cr
#'   \code{measure_ets_vcov} generates the variance-covariance matrix for 
#'   the regressors.
#'
#' @param etas Matrix of regressors (columns) across iterations of the 
#'   sampler (rows), as returned from \code{\link{est_regressors}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{summarize_etas}: table of summary statistics for chunk-level
#'   regressors including mean, median, mode, posterior interval, standard
#'   deviation, MCMC error, autocorrelation, and effective sample size for 
#'   each regressor. \cr \cr
#'   \code{measure_eta_vcov}: variance-covariance matrix for chunk-level
#'   regressors.
#'
#' @export 
#'
summarize_etas <- function(etas, control = list()){
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

#'
#' @rdname summarize_etas
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
#'   the sampler (rows) or \code{NULL} if no change points are in the model,
#'   as returned from \code{\link{est_change points}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   Time Series model. Values not input assume defaults set by 
#'   \code{\link{TS_control}}.
#'
#' @return \code{summarize_rhos}: table of summary statistics for change point
#'   locations including mean, median, mode, posterior interval, standard
#'   deviation, MCMC error, autocorrelation, and effective sample size for 
#'   each change point location. \cr \cr
#'   \code{measure_rho_vcov}: variance-covariance matrix for change 
#'   point locations.
#' 
#' @export 
#'
summarize_rhos <- function(rhos, control = list()){
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

#'
#' @rdname summarize_rhos
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
