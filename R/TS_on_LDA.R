
TS_on_LDA <- function(LDA_models, document_covariate_table, formulas = ~ 1, 
                      nchangepoints = 0, timename = "time", weights = NULL, 
                      control = list()){
  check_TS_on_LDA_inputs(LDA_models, document_covariate_table, formulas, 
                         nchangepoints, timename, weights, control)
  control <- do.call("TS_control", control)
  mods <- expand_TS(LDA_models, formulas, nchangepoints)
  nmods <- nrow(mods)
  TSmods <- vector("list", nmods)

  for(i in 1:nmods){
    print_model_run_message(mods, i, LDA_models, control)
    formula_i <- mods$formula[[i]]
    nchangepoints_i <- mods$nchangepoints[i]
    data_i <- prep_TS_data(document_covariate_table, LDA_models, mods, i)
    TSmods[[i]] <- TS(data_i, formula_i, nchangepoints_i, timename, weights, 
                      control)
  }
  package_TS_on_LDA(TSmods, LDA_models, mods)

}
prep_TS_data <- function(document_covariate_table, LDA_models, mods, i = 1){
  check_document_covariate_table(document_covariate_table, LDA_models)
  check_LDA_models(LDA_models)
  if(is(LDA_models, "LDA")){
    LDA_models <- c(LDA_models)
    class(LDA_models) <- c("LDA_set", "list")
  }
  data_i <- document_covariate_table
  data_i$gamma <- LDA_models[[mods$LDA[i]]]@gamma
  data_i
}

select_TS <- function(TS_models, control = list()){
  if (!("TS_on_LDA" %in% class(TS_models))){
    stop("TS_models must be of class TS_on_LDA")
  }
  check_control(control)
  control <- do.call("TS_control", control)
  measurer <- control$measurer
  selector <- control$selector
  TS_measured <- vapply(TS_models, measurer, 0) %>%
                  matrix(ncol = 1)
  TS_selected <- apply(TS_measured, 2, selector) 
  which_selected <- which(TS_measured %in% TS_selected)
  if (length(which_selected) > 1){
    warning("Selection results in multiple models, returning first")
    which_selected <- which_selected[1]
  }
  out <- TS_models[[which_selected]]
  class(out)  <- c("TS_fit", "list") 
  out
}


package_TS_on_LDA <- function(TSmods, LDA_models, models){
  check_LDA_models(LDA_models)
  if(is(LDA_models, "LDA")){
    LDA_models <- c(LDA_models)
    class(LDA_models) <- c("LDA_set", "list")
  }
  nmodels <- nrow(models)
  nms <- rep(NA, nmodels)
  for (i in 1:nmodels){
    nms[i] <- paste0(names(LDA_models)[models$LDA[i]], ", ", 
                     deparse(models$formula[[i]]), ", ", 
                     models$nchangepoints[i], " changepoints")
  }
  names(TSmods) <- nms
  class(TSmods) <- list("TS_on_LDA", "list")
  TSmods
}


print.TS_on_LDA <- function(x, ...){
  print(names(x))
}


print_model_run_message <- function(models, i, LDA_models, control){
  control <- do.call("TS_control", control)
  equation <- deparse(models$formula[[i]])
  chngpt_msg <- paste0("with ", models$nchangepoints[i], " changepoints ")
  reg_msg <- paste0("and equation ", equation)
  ts_msg <- paste0(chngpt_msg, reg_msg)
  lda_msg <- names(LDA_models)[models$LDA[i]]
  msg <- paste0("Running TS model ", ts_msg, " on LDA model ", lda_msg, "\n")
  messageq(msg, control$quiet)
}

expand_TS <- function(LDA_models, formulas, nchangepoints){
  check_LDA_models(LDA_models)
  check_nchangepoints(nchangepoints)
  if (is(LDA_models, "LDA")) {
    LDA_models <- c(LDA_models)
    class(LDA_models) <- c("LDA_set", "list")
  }
  if (!is(formulas, "list")) {
    if (is(formulas, "formula")) {
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain formula(s)")
    }
  } else if (!all(vapply(formulas, is, TRUE, "formula"))) {
      stop("formulas does not contain all formula(s)")
  }
  formulas
  
  out <- formulas
  for (i in seq_along(formulas)) {
    tformula <- paste(as.character(formulas[[i]]), collapse = "")
    out[[i]] <- as.formula(paste("gamma", tformula))
  }
  formulas <- out
  nmods <- length(LDA_models)
  mods <- 1:nmods
  out <- expand.grid(mods, formulas, nchangepoints, stringsAsFactors = FALSE)
  colnames(out) <- c("LDA", "formula", "nchangepoints") 
  out
}

check_nchangepoints <- function(nchangepoints){
  if (!is.numeric(nchangepoints) || any(nchangepoints %% 1 != 0)){
    stop("nchangepoints must be integer-valued")
  }
  if (any(nchangepoints < 0)){
    stop("nchangepoints must be non-negative")
  }
  return()
}

check_weights <- function(weights){
  if(is.logical(weights)){
    if(weights){
      return()
    } else{
      stop("if logical, weights need to be TRUE")
    }   
  }
  if(!is.null(weights)){
    if (!is.numeric(weights)){
      stop("weights vector must be numeric")
    }
    if (any(weights <= 0)){
      stop("weights must be positive")
    }
    if (round(mean(weights)) != 1){
      warning("weights should have a mean of 1, fit may be unstable")
    }
  }
  return()
}


check_LDA_models <- function(LDA_models){
  if(("LDA_set" %in% class(LDA_models)) == FALSE){
    if(is(LDA_models, "LDA") == FALSE){
      stop("LDA_models is not an LDA object or LDA_set object")
    }
  }
  return()
}

check_document_covariate_table <- function(document_covariate_table, 
                                           LDA_models = NULL,
                                           document_term_table = NULL){
  dct_df <- tryCatch(data.frame(document_covariate_table),
                     warning = function(x){NA}, error = function(x){NA})
  if(is(LDA_models, "LDA")){
    LDA_models <- c(LDA_models)
    class(LDA_models) <- c("LDA_set", "list")
  }
  if (length(dct_df) == 1 && is.na(dct_df)){
    stop("document_covariate_table is not conformable to a data frame")
  }
  if (!is.null(LDA_models)){
    if (nrow(data.frame(document_covariate_table)) != 
        nrow(LDA_models[[1]]@gamma)){
      stop("number of documents in covariate table is not equal to number of 
        documents observed")
    }
  } else if (!is.null(document_term_table)){
    if (nrow(data.frame(document_covariate_table)) != 
        nrow(data.frame(document_term_table))){
      stop("number of documents in covariate table is not equal to number of 
        documents observed")
    }
  }
  return()
}

check_timename <- function(document_covariate_table, timename){
  if (!("character" %in% class(timename))){
    stop("timename is not a character value")
  }
  if (length(timename) > 1){
    stop("timename can only be one value")
  }
  covariate_names <- colnames(document_covariate_table)
  if ((timename %in% covariate_names) == FALSE){
    stop("timename not present in document covariate table")
  }
  time_covariate <- document_covariate_table[ , timename]
  if (!(is.Date(time_covariate)) & 
      (!is.numeric(time_covariate) || !all(time_covariate %% 1 == 0))){
    stop("covariate indicated by timename is not an integer or a date")
  }
  return()
}


check_formulas <- function(formulas, document_covariate_table, 
                           control = list()){
  check_document_covariate_table(document_covariate_table)
  check_control(control)
  control <- do.call("TS_control", control)
  # response <- control$response
  dct <- document_covariate_table
  if (!is(formulas, "list")) {
    if (is(formulas, "formula")) {
      formulas <- c(formulas)
    } else{
      stop("formulas does not contain formula(s)")
    }
  } else if (!all(vapply(formulas, is, TRUE, "formula"))) {
      stop("formulas does not contain all formula(s)")
  }
  resp <- unlist(lapply(lapply(formulas, terms), attr, "response"))
  pred <- unlist(lapply(lapply(formulas, terms), attr, "term.labels"))
  if (any(resp != 0)) {
    stop("formula inputs should not include response variable")
  }
  if (!all(pred %in% colnames(dct))) {
    misses <- pred[which(pred %in% colnames(dct) == FALSE)]
    mis <- paste(misses, collapse = ", ")
    stop(paste0("formulas include predictors not present in data: ", mis))
  }
  return()
}



check_TS_on_LDA_inputs <- function(LDA_models, document_covariate_table, 
                            formulas = ~ 1, nchangepoints = 0,  
                            timename = "time", weights = NULL,
                            control = list()){
  check_LDA_models(LDA_models)
  check_document_covariate_table(document_covariate_table, LDA_models)
  check_timename(document_covariate_table, timename)
  check_formulas(formulas, document_covariate_table, control)  
  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_control(control)
}
