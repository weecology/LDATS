meta_LDA_TS <- function(data, topics = 2, nseeds = 1, formulas = ~ 1,
                        nchangepoints = 0, niterations = 1, timename = "time",  
                        weights = TRUE, control = list()){

  check_meta_LDA_TS_inputs(data, topics, nseeds, formulas, nchangepoints,  
                      timename, weights, control)
  control <- do.call("meta_LDA_TS_control", control)
  data <- conform_LDA_TS_data(data, control$quiet)

  mLDATS <- vector("list", length = niterations)
  for(i in 1:niterations){

# data selection will happen here

    LDATS_fit <- LDA_TS(data = data, topics = topics, nseeds = nseeds,
                        formulas = formulas, nchangepoints = nchangepoints,
                        timename = timename, weights = weights, 
                        control = control)
    mLDATS[[i]] <- list(LDA_TS = LDATS_fit, data = data)
  }
  names(mLDATS) <- paste0("meta_LDA_TS iteration ", 1:niterations)



# we're going to want to measure, select (although not using) and package
  mLDATS
}



check_meta_LDA_TS_inputs <- function(data = NULL,
                                     topics = 2, nseeds = 1, formulas = ~ 1, 
                                     nchangepoints = 0,
                                     timename = "time", 
                                     weights = TRUE, 
                                     control = list()){
  check_control(control)
  control <- do.call("meta_LDA_TS_control", control)
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


meta_LDA_TS_control <- function(quiet = FALSE, measurer_LDA = AIC, 
                                selector_LDA = min, iseed = 2,
                                memoise = TRUE, response = "gamma", 
                                lambda = 0, measurer_TS = AIC, 
                                selector_TS = min, ntemps = 6, 
                                penultimate_temp = 2^6, ultimate_temp = 1e10, 
                                q = 0, nit = 1e4, magnitude = 12, burnin = 0, 
                                thin_frac = 1, summary_prob = 0.95, 
                                seed = NULL, measurer_LDA_TS = logLik, 
                                selector_LDA_TS = NULL, ...){

  list(quiet = quiet, measurer_LDA = measurer_LDA, 
       selector_LDA = selector_LDA,
       iseed = iseed, memoise = memoise, response = response, lambda = lambda, 
       measurer_TS = measurer_TS, selector_TS = selector_TS, ntemps = ntemps, 
       penultimate_temp = penultimate_temp, ultimate_temp = ultimate_temp, 
       q = q, nit = nit, magnitude = magnitude, burnin = burnin, 
       thin_frac = thin_frac, summary_prob = summary_prob, seed = seed,
       measurer_LDA_TS = measurer_LDA_TS, selector_LDA_TS = selector_LDA_TS)
}