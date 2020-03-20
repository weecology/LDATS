##
#
#  working in here 
#
#  looks pretty nice!
#
#  generally pretty well working, needs to get tidied up considerably tho
#  pull redundant code out into functions etc
#   and figure out better namings etc
#  also make sure when no changepoint that it still says estimating regressor



TS <- function(LDAs, data, formulas = ~ 1, nchangepoints = 0, 
               timename = "time", weights = NULL, control = list()){
  control <- do.call("TS_control", control)
  messageq("----- Time Series Analyses -----", control$quiet)
  TSs <- prep_TS_models(LDAs = LDAs, data = data, formulas = formulas,
                        nchangepoints = nchangepoints, timename = timename,
                        weights = weights, control = control)
  nTS <- length(TSs)
  for (i in 1:nTS){
    TSs[[i]] <- sequential_TS(TS = TSs[[i]], control = control)
  }
  selected_TSs <- select_TS(TSs = TSs, control = control)
  package_TS(selected_TSs = selected_TSs, TSs = TSs, control = control)
}



select_TS <- function(TSs, control = list()){

  vals <- measure_TS(TSs = TSs, control = control)
  fun <- control$selector_function
  args <- update_list(control$selector_args, x = vals)
  selection <- do.call(what = fun, args = args)
  TSs[selection]  
}

measure_TS <- function(TSs, control = list()){
  fun <- control$measurer_function
  args <- control$measurer_args
  nTSs <- length(TSs)
  vals <- rep(NA, nTSs)
  for(i in 1:nTSs){
    args <- update_list(args, object = TSs[[i]])
    vals_i <- do.call(what = fun, args = args)
    if(length(vals_i) != 0){
      vals[i] <- vals_i
    }
  }
  vals
}



package_TS <- function(selected_TSs, TSs, control = list()){
  out <- list(selected_TSs = selected_TSs, TSs = TSs, control = control)
  class(out) <- c("TS_set", "list")
  out
}






TS_control <- function(response = "multinom",
                       method = "ldats_classic",
                       method_args = list(ntemps = 6, penultimate_temp = 2^6, 
                                          ultimate_temp = 1e10, q = 0, 
                                          nit = 1e4, magnitude = 12, 
                                          burnin = 0, thin_frac = 1, 
                                          memoise = TRUE,
                                          quiet = FALSE),
                       summary_prob = 0.95,
                        measurer_function = AIC,
                        measurer_args = list(),
                        selector_function = which.min,
                        selector_args = list(), 
                       soften = TRUE, 
                       quiet = FALSE){
  list(response = response, method = method, method_args = method_args, 
       measurer_function = measurer_function, measurer_args = measurer_args, 
       selector_function = selector_function, selector_args = selector_args,
       summary_prob = summary_prob, soften = soften, quiet = quiet)
}


sequential_TS <- function(TS, control = list()){
  control <- do.call("TS_control", control)
  TS_msg(TS = TS, quiet = control$quiet)
  rho_dist <- est_changepoints(TS = TS, control = control)
  eta_dist <- est_regressors(rho_dist = rho_dist, TS = TS, control = control)

  package_sequential_TS(TS = TS, rho_dist = rho_dist, eta_dist = eta_dist, 
                        control = control)
}

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






print.TS <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !(names(x) %in% hid)
  print(x[notHid])
}


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

measure_eta_vcov <- function(etas){
  if (!is.matrix(etas)){
    stop("expecting etas to be a matrix")
  }
  out <- var(etas)
  colnames(out) <- colnames(etas)
  rownames(out) <- colnames(etas)
  out
}
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




est_changepoints <- function(TS, control = list()){
  if (TS$nchangepoints == 0){
    return(NULL)
  }
  fun <- control$method
  args <- list(TS = TS, control = control$method_args)
  soft_call(fun = fun, args = args, soften = control$soften)
}


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


prep_TS_models <- function(LDAs, data, formulas = ~ 1, nchangepoints = 0, 
                           timename = "time", weights = NULL, 
                           control = list()){

  if (!is(formulas, "list")) {
    if (is(formulas, "formula")) {
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain formula(s)")
    }
  } else if (!all(vapply(formulas, is, TRUE, "formula"))) {
      stop("formulas does not contain all formula(s)")
  }

  formulas2 <- formulas
  for (i in seq_along(formulas)) {
    tformula <- paste(as.character(formulas[[i]]), collapse = "")
    formulas2[[i]] <- as.formula(paste("gamma", tformula))
  }
  formulas <- formulas2
  nmods <- length(LDAs[[1]])
  mods <- 1:nmods
  tab <- expand.grid(mods, formulas, nchangepoints, stringsAsFactors = FALSE)
  colnames(tab) <- c("LDA", "formula", "nchangepoints") 
  
  nTSs <- NROW(tab)
  TSs <- vector("list", length = nTSs)
  for(i in 1:nTSs){
    lda <- LDAs[[1]][[tab$LDA[i]]]
 
    ts_data <- lda$data
    ts_data$train$ts_data <- lda$data$train$document_covariate_table
    ts_data$train$ts_data$gamma <- lda$document_topic_matrix

    ts_data$test$ts_data <- lda$data$test$document_covariate_table
    ts_data$test$ts_data$gamma <- lda$test_document_topic_matrix

    weights <- iftrue(weights, 
                      document_weights(lda$data$train$document_term_table))

    TSs[[i]] <- list(data = ts_data,
                     data_subset = lda[["data_subset"]],
                     formula = tab$formula[[i]],
                     nchangepoints = tab$nchangepoints[i], 
                     weights = weights,
                     timename = timename,
                     response = control$response,
                     topics = lda$topics, rep = lda$rep)
  }
  name_tab <- data.frame(paste("LDA", tab[ , 1]), 
                         paste(",", tab[ , 2]),
                         paste(",", tab[ , 3], "changepoints"))
  names(TSs) <- apply(name_tab, 1, paste0, collapse = "")
  TSs
}



TS_msg <- function(TS, quiet = FALSE){
  subset_msg <- paste0("  - data subset ", TS$data_subset)
  topic_msg <- paste0(", ", TS$topics, " topics")
  rep_msg <- paste0(", replicate ", TS$rep)

  formula_msg <- paste0(", ", deparse(TS$formula))
  nchangepoints <- TS$nchangepoints
  txt <- ifelse(nchangepoints == 1, " changepoint", " changepoints")
  changepoints_msg <- paste0(", ", nchangepoints, txt)
  msg <- paste0(subset_msg, topic_msg, rep_msg, formula_msg, changepoints_msg)
  messageq(msg, quiet)
}



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
  return()
}


prep_pbar <- function(control = list(), bar_type = "rho", nr = NULL){
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
    msg <- "  Estimating changepoint distribution"
    out <- progress_bar$new(form, control$nit, width = 60)
  }
  if (bar_type == "eta"){
    msg <- "  Estimating regressor distribution"
    out <- progress_bar$new(form, nr, width = 60)
  }
  messageq(msg, control$quiet)
  out
}

update_pbar <- function(pbar, control = list()){
  if (!("progress_bar" %in% class(pbar))){
    stop("pbar must be of class progress_bar")
  }
  if (control$quiet){
    return()
  }
  pbar$tick()
}

logLik.TS <- function(object, ...){
  val <- object$logLik
  attr(val, "df") <- object$nparams
  attr(val, "nobs") <- nrow(object$data)
  class(val) <- "logLik"
  val
}


#' @title Check that a set of change point locations is proper
#' 
#' @description Check that the change point locations are \code{numeric}
#'   and conformable to \code{interger} values. 
#'   
#' @param changepoints Change point locations to evaluate.
#' 
#' @return An error message is thrown if \code{changepoints} are not proper,
#'   else \code{NULL}.
#'
#' @examples
#'   check_changepoints(100)
#'
#' @export
#'
check_changepoints <- function(changepoints = NULL){
  if (is.null(changepoints)){
    return()
  }
  if (!is.numeric(changepoints) || any(changepoints %% 1 != 0)){
    stop("changepoints must be integer-valued")
  }
}


#' @title Prepare the time chunk table for a multinomial change point 
#'   Time Series model
#'
#' @description Creates the table containing the start and end times for each
#'   chunk within a time series, based on the change points (used to break up
#'   the time series) and the range of the time series. If there are no 
#'   change points (i.e. \code{changepoints} is \code{NULL}, there is still a
#'   single chunk defined by the start and end of the time series.
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables, but specifically here containing the column indicated
#'   by the \code{timename} input. 
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points. Must be conformable to \code{integer} values. 
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @return \code{data.frame} of \code{start} and \code{end} times (columns)
#'   for each chunk (rows).
#'
#' @examples
#'   data(rodents)
#'   dtt <- rodents$document_term_table
#'   lda <- LDA_set(dtt, 2, 1, list(quiet = TRUE))
#'   dct <- rodents$document_covariate_table
#'   dct$gamma <- lda[[1]]@gamma
#'   chunks <- prep_chunks(dct, changepoints = 100, timename = "newmoon")   
#'
#' @export 
#'
prep_chunks <- function(data, changepoints = NULL, 
                        timename = "time"){
  start <- c(min(data[ , timename]), changepoints + 1)   
  end <- c(changepoints, max(data[ , timename])) 
  data.frame(start, end)
}


#' @title Verify the change points of a multinomial time series model
#'
#' @description Verify that a time series can be broken into a set 
#'   of chunks based on input change points. 
#'
#' @param data Class \code{data.frame} object including the predictor and 
#'   response variables.
#'
#' @param changepoints Numeric vector indicating locations of the change 
#'   points. Must be conformable to \code{integer} values. 
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior.
#'
#' @return Logical indicator of the check passing \code{TRUE} or failing
#'   \code{FALSE}.
#'
#' @examples
#'   data(rodents)
#'   dtt <- rodents$document_term_table
#'   lda <- LDA_set(dtt, 2, 1, list(quiet = TRUE))
#'   dct <- rodents$document_covariate_table
#'   dct$gamma <- lda[[1]]@gamma
#'   verify_changepoint_locations(dct, changepoints = 100, 
#'                                timename = "newmoon")   
#'
#' @export 
#'
verify_changepoint_locations <- function(data, changepoints = NULL, 
                                     timename = "time"){

  if (is.null(changepoints)){
    return(TRUE)
  }

  first_time <- min(data[ , timename])
  last_time <- max(data[ , timename])
  time_check <- any(changepoints <= first_time | changepoints >= last_time)
  sort_check <- is.unsorted(changepoints, strictly = TRUE)

  !(time_check | sort_check)
}
