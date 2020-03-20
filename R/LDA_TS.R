# pulled back to just what's working so far
# the goal is to have a single LDA_TS function, no need for multi etc
#

LDA_TS <- function(data, topics = 2, reps = 1, formulas = ~ 1,
                   nchangepoints = 0, timename = "time", weights = TRUE, 
                   control = list()){
  control <- do.call("LDA_TS_control", control)
  data <- conform_data(data = data, control = control)
  LDAs <- LDA(data = data, topics = topics, reps = reps, 
              control = control$LDA_control)
  TSs <- TS(LDAs = LDAs, data = data, formulas = formulas, 
            nchangepoints = nchangepoints, timename = timename, 
            weights = weights, control = control$TS_control) 
 
}



LDA_TS_control <- function(LDA_model = topicmodels_LDA, 
                           LDA_model_args = 
                             list(method = "VEM", seeded = TRUE),
                           LDA_measurer = AIC,
                           LDA_measurer_args = list(),
                           LDA_selector = which.min,
                           LDA_selector_args = list(), 
                           TS_model = sequential_TS,
                           TS_model_args = list(),
                           TS_response = "multinom",
                           TS_response_args = list(),
                           TS_method = "ldats_classic",
                           TS_method_args = ldats_classic_control(),
                           TS_measurer = AIC,
                           TS_measurer_args = list(),
                           TS_selector = which.min,
                           TS_selector_args = list(),
                           nsubsets = 1,
                           subset_rule = NULL,
                           summary_prob = 0.95,
                           soften = TRUE, 
                           quiet = FALSE, ...){

  LDA_control <- LDA_control(model = LDA_model, model_args = LDA_model_args, 
                             measurer = LDA_measurer, 
                             measurer_args = LDA_measurer_args, 
                             selector = LDA_selector, 
                             selector_args = LDA_selector_args,
                             nsubsets = nsubsets, subset_rule = subset_rule,
                             soften = soften, quiet = quiet)

  TS_control <- TS_control(response = TS_response, 
                           response_args = TS_response_args, 
                           method = TS_method, 
                           method_args = TS_method_args, 
                           measurer = TS_measurer,
                           measurer_args = TS_measurer_args, 
                           selector = TS_selector, 
                           selector_args = TS_selector_args,
                           summary_prob = summary_prob, 
                           soften = soften, quiet = quiet)


  list(LDA_control = LDA_control, TS_control = TS_control, 
       nsubsets = nsubsets, subset_rule = subset_rule,
       soften = soften, quiet = quiet)
}


