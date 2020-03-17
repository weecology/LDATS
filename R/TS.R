##
#
#  working in here 




TS <- function(LDAs, data, formulas = ~ 1, nchangepoints = 0, 
               timename = "time", weights = NULL, control = list()){
  control <- do.call("TS_control", control)
  messageq("----Time Series Analyses----", control$quiet)
  TSs <- prep_TS_models(LDAs = LDAs, data = data, formulas = formulas,
                        nchangepoints = nchangepoints, timename = timename,
                        weights = weights, control = control)
  nTS <- length(TSs)
  for (i in 1:nTS){
    TSs[[i]] <- TS_call(TS = TSs[[i]], control = control)
  }
  selected_TSs <- select_TS(TSs = TSs, control = control)
  package_TS(selected_TSs = selected_TSs, TSs = TSs, control = control)
}




# the defining characteristic here is the sequential
#  estimation of rho then eta, ostensibly it could be done together


sequential_TS <- function(TS, control = list()){
  rho_dist <- est_changepoints(TS = TS, control = control)
  eta_dist <- est_regressors(rho_dist = rho_dist, TS = TS, control = control)

}

est_changepoints <- function(TS, control = list()){

  if (nchanTS$nchangepoints == 0){
    return(NULL)
  }

x <- memoizer(TS, control = control$TS_fit_args)

   
}




memoizer <- function(TS, control = list()){
    
  data <- TS$data$train$ts_data
  nchangepoints <- TS$nchangepoints
  formula <- TS$formula
  weights <- TS$weights
  timename <- TS$timename


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

    TSs[[i]] <- list(data = ts_data,
                     data_subset = lda[["data_subset"]],
                     formula = tab$formula[[i]],
                     nchangepoints = tab$nchangepoints[i], 
                     weights = weights,
                     timename = timename)
  }
  name_tab <- data.frame(paste("LDA", tab[ , 1]), 
                         paste(",", tab[ , 2]),
                         paste(",", tab[ , 3], "changepoints"))
  names(TSs) <- apply(name_tab, 1, paste0, collapse = "")
  TSs
}



TS_control <- function(TS_function = multinom_TS, 
                       TS_args = list(),
                       TS_fit_function = memoizer,
                       TS_fit_args = list(ntemps = 6, penultimate_temp = 2^6, 
                                          ultimate_temp = 1e10, q = 0, 
                                          nit = 1e4, magnitude = 12, 
                                          burnin = 0, thin_frac = 1, 
                                          quiet = FALSE),
                       soften = TRUE, 
                       quiet = FALSE){
  list(TS_function = TS_function, 
       TS_args = TS_args, 
       TS_fit_function = TS_fit_function, TS_fit_args = TS_fit_args, 
       soften = soften, quiet = quiet)
}




TS_call <- function(TS = NULL, control = list()){
  control <- do.call("TS_control", control)  
  TS_msg(TS = TS, quiet = control$quiet)
  fun <- control$TS_function
  args <- update_list(control$TS_args, TS = TS)
  if(control$soften){
    tryCatch(do.call(what = fun, args = args), 
             warning = function(x){eval(x$call)}, 
             error = function(x = list()){list(error = x$message)})
  } else{
    do.call(what = fun, args = args)
  }
}


TS_msg <- function(TS, quiet = FALSE){
  messageq("hi how are you?", quiet)
}















predict.TS_fit <- function(object, newdata = NULL, control = list(), ...){
  if(is.null(newdata)){
    newdata <- object$data
  }
  control <- do.call("TS_control", control)
  nit <- object$control$nit
  rhos <- object$rhos
  etas <- object$etas
  nnewdata <- NROW(newdata)
  formula <- object$formula
  out <- matrix(NA, nrow = nit, ncol = nnewdata)
  for(i in 1:nit){
    
    out[i , ] <- predicts
  }

}

TSx <- function(data, formula = gamma ~ 1, nchangepoints = 0, 
               timename = "time", weights = NULL, control = list()){
  check_TS_inputs(data, formula, nchangepoints, timename, weights, control)
  control <- do.call("TS_control", control)
  set.seed(control$seed)
  data <- time_order_data(data, timename = timename)
  rho_dist <- est_changepoints(data, formula, nchangepoints, timename, 
                               weights, control)
  eta_dist <- est_regressors(rho_dist, data, formula, timename, weights, 
                             control)
  package_TS(data, formula, timename, weights, control, rho_dist, eta_dist)
}

check_TS_inputs <- function(data, formula = gamma ~ 1, nchangepoints = 0, 
                            timename = "time", weights = NULL, 
                            control = list()){
  check_formula(data, formula)  
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_timename(data, timename)
  check_control(control)
  return() 
}

package_TSx <- function(data, formula, timename, weights, control, rho_dist, 
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
  AIC <- -2 * logLik + 2 * nparams

  out <- list(data = data, formula = formula, nchangepoints = nchangepoints,
              timename = timename, weights = weights,
              control = control, lls = lls, rhos = rhos,
              etas = eta_dist, ptMCMC_diagnostics = ptMCMC_diagnostics,
              rho_summary = rho_summary, rho_vcov = rho_vcov,
              eta_summary = eta_summary, eta_vcov = eta_vcov,
              logLik = logLik, nparams = nparams, AIC = AIC)
  class(out) <- c("TS_fit", "list")
  to_hide <- c("data", "weights", "control", "lls", "rhos", "etas", 
               "rho_vcov", "eta_vcov")
  if (nchangepoints == 0){
    to_hide <- c(to_hide, "ptMCMC_diagnostics", "rho_summary")
  }
  attr(out, "hidden") <- to_hide
  out
}
print.TS_fit <- function(x, ...){
  hid <- attr(x, "hidden")
  notHid <- !names(x) %in% hid
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

est_regressorsx <- function(rho_dist, data, formula, timename, weights, 
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


est_changepointsx <- function(data, formula, nchangepoints, timename, weights, 
                             control = list()){
  check_TS_inputs(data, formula, nchangepoints, timename, weights, control)
  control <- do.call("TS_control", control)
  data <- time_order_data(data, timename = timename)
  if (nchangepoints == 0){
    return(NULL)
  }
  saves <- prep_saves(nchangepoints, control)

# break this into a "classic" approach or whatever and give it its own
#  function that is akin to temper

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
# the separate function runs to here
# and then process the output

  process_saves(saves, control)
}

prep_pbar <- function(control = list(), bar_type = "rho", 
                      nr = NULL){
  check_control(control)
  control <- do.call("TS_control", control)
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
  check_control(control)
  control <- do.call("TS_control", control)
  if (control$quiet){
    return()
  }
  pbar$tick()
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

TS_controlx <- function(memoise = TRUE, response = "gamma", lambda = 0, 
                       measurer = AIC, selector = min, ntemps = 6, 
                       penultimate_temp = 2^6, ultimate_temp = 1e10, q = 0, 
                       nit = 1e4, magnitude = 12, quiet = FALSE, burnin = 0, 
                       thin_frac = 1, summary_prob = 0.95, seed = NULL,
                       model_fun = "multinom"){
  list(memoise = memoise, response = response, lambda = lambda, 
       measurer = measurer, selector = selector, ntemps = ntemps, 
       penultimate_temp = penultimate_temp, ultimate_temp = ultimate_temp, 
       q = q, nit = nit, magnitude = magnitude, quiet = quiet, 
       burnin = burnin, thin_frac = thin_frac, summary_prob = summary_prob,
       seed = seed, model_fun = model_fun)

}

logLik.TS_fit <- function(object, ...){
  val <- object$logLik
  attr(val, "df") <- object$nparams
  attr(val, "nobs") <- nrow(object$data)
  class(val) <- "logLik"
  val
}