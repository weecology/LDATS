# runs multiple iterations of an LDA_TS model (or model suite) where each time
# a subset of the data are held out

multi_LDA_TS <- function(data, topics = 2, nseeds = 1, formulas = ~ 1,
                         nchangepoints = 0, niterations = 1, 
                         timename = "time",  
                         weights = TRUE, control = list()){

  check_multi_LDA_TS_inputs(data, topics, nseeds, formulas, nchangepoints,  
                            timename, weights, control)
  control <- do.call("multi_LDA_TS_control", control)
  data <- conform_LDA_TS_data(data, control$quiet)

  mLDATS <- vector("list", length = niterations)
  for(i in 1:niterations){

    data_i <- split_LDA_TS_data(data, control, i)   
    LDATS_fit <- LDA_TS(data = data_i[["train"]], topics = topics, 
                        nseeds = nseeds,
                        formulas = formulas, nchangepoints = nchangepoints,
                        timename = timename, weights = weights, 
                        control = control)
    test <- test_LDA_TS(LDATS_fit, data = data_i[["test"]], control)
    mLDATS[[i]] <- list(LDA_TS = LDATS_fit, test = test, data = data_i)
  }
  names(mLDATS) <- paste0("multi_LDA_TS iteration ", 1:niterations)
  
  # select and package
  
  mLDATS
}

test_LDA_TS <- function(LDATS_fit, data, control = list()){
  # need to make predictions cleanly
  # predict method on LDA_TS which combines predict method on TS s
  pred_data <- predict(LDATS_fit, newdata = data, control = control) 
  vapply(TS_models, measurer, 0)
}

# rule should be a function, yeah?
#  yes
# if rule is NULL, nothing is held out, like with a classical ldats
# its really test train out

split_LDA_TS_data <- function(data = NULL, control = list(), iteration = 1){
  control <- do.call("multi_LDA_TS_control", control)
  data <- conform_LDA_TS_data(data, control$quiet)
  dtt <- data[["document_term_table"]]
  dct <- data[["document_covariate_table"]]
  rule <- control$rule
  if(is.null(rule)){
    rule <- null_rule
  }
  arglist <- list(data = dtt, iteration = iteration)
  test_train <- do.call(rule, arglist)

  train <- list(document_term_table = dtt[test_train == "train", ],
                document_covariate_table = dct[test_train == "train", ])
  test <- list(document_term_table = dtt[test_train == "test", ],
               document_covariate_table = dct[test_train == "test", ])
  list(train = train, test = test, test_train = test_train)
}


null_rule <- function(data, iteration = 1){
  n <- NROW(data)
  rep("train", n)
}

# simple leave one outs with no buffer

#  for use with exhaustive approach,
#  as it assumes 1:1 between iteraction and datum location to drop

exhaust_loo <- function(data, iteration = 1){
  leave_p_out(data = data, random = FALSE, locations = iteration)
}

random_loo <- function(data, iteration = 1){
  leave_p_out(data = data)
}

# fully flexible leave p out function allowing for buffers
# if random the test data are selected randomly, otherwise locations are used
leave_p_out <- function(data, p = 1, pre = 0, post = 0, 
                        random = TRUE, locations = NULL){
  n <- NROW(data)
  test_train <- rep("train", n)

  if(random){
    locations <- sample(1:n, p)
  }

  for(i in 1:p){
    hold_out <- (locations[i] - pre):(locations[i] + post)
    test_train[hold_out] <- "out"
  }
  test_train[locations] <- "test"
  test_train
}




check_multi_LDA_TS_inputs <- function(data = NULL,
                                     topics = 2, nseeds = 1, formulas = ~ 1, 
                                     nchangepoints = 0,
                                     timename = "time", 
                                     weights = TRUE, 
                                     control = list()){
  check_control(control)
  control <- do.call("multi_LDA_TS_control", control)
  data <- conform_LDA_TS_data(data)
  weights <- iftrue(weights, document_weights(data$document_term_table))
  check_document_covariate_table(data$document_covariate_table, 
                               document_term_table = data$document_term_table)
  check_timename(data$document_covariate_table, timename)
  ts_control <- TS_control(control)
  check_formulas(formulas, data$document_covariate_table, ts_control)  

  check_nchangepoints(nchangepoints)
  check_weights(weights)
  check_document_term_table(data$document_term_table)
  check_topics(topics)
  check_seeds(nseeds)

}

# note that the controls here override the standard LDA_TS inputs with NULL
#  for each of the measurer and selector functions for LDA and TS
#  default here is I for selector, which just passes along (identity function)
#
multi_LDA_TS_control <- function(quiet = FALSE, measurer_LDA = logLik, 
                                selector_LDA = I, iseed = 2,
                                memoise = TRUE, response = "gamma", 
                                lambda = 0, measurer_TS = logLik, 
                                selector_TS = I, ntemps = 6, 
                                penultimate_temp = 2^6, ultimate_temp = 1e10, 
                                q = 0, nit = 1e4, magnitude = 12, burnin = 0, 
                                thin_frac = 1, summary_prob = 0.95, 
                                seed = NULL, measurer_LDA_TS = logLik, 
                                selector_LDA_TS = I, 
                                rule = exhaust_loo, ...){

  list(quiet = quiet, measurer_LDA = measurer_LDA, 
       selector_LDA = selector_LDA,
       iseed = iseed, memoise = memoise, response = response, lambda = lambda, 
       measurer_TS = measurer_TS, selector_TS = selector_TS, ntemps = ntemps, 
       penultimate_temp = penultimate_temp, ultimate_temp = ultimate_temp, 
       q = q, nit = nit, magnitude = magnitude, burnin = burnin, 
       thin_frac = thin_frac, summary_prob = summary_prob, seed = seed,
       measurer_LDA_TS = measurer_LDA_TS, selector_LDA_TS = selector_LDA_TS,
       rule = rule)
}