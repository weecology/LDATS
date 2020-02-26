# pulled back to just what's working so far
# the goal is to have a single LDA_TS function, no need for multi etc
#

LDA_TS <- function(data, topics = 2, reps = 1, formulas = ~ 1,
                   nchangepoints = 0, timename = "time", weights = TRUE, 
                   control = list()){
  control <- do.call("LDA_TS_control", control)
  data <- conform_data(data = data, control = control)
  LDAs <- LDA(data = data, topics = topics, nseeds = nseeds, 
              control = control$LDA_control)
  

}



LDA_TS_control <- function(LDA_function = topicmodels_LDA, 
                        LDA_args = list(method = "VEM", seeded = TRUE),
                        LDA_measurer_function = AIC,
                        LDA_measurer_args = list(),
                        LDA_selector_function = which.min,
                        LDA_selector_args = list(), 
                        nsubsets = 1,
                        subset_rule = NULL,
                        soften = TRUE, 
                        quiet = FALSE){

  LDA_control <- LDA_control(LDA_function = LDA_function, LDA_args = LDA_args, 
                            measurer_function = LDA_measurer_function, 
                            measurer_args = LDA_measurer_args, 
                            selector_function = LDA_selector_function, 
                            selector_args = LDA_selector_args,
                            nsubsets = nsubsets, subset_rule = subset_rule,
                            soften = soften, quiet = quiet)
  TS_control <- TS_control()
  list(LDA_control = LDA_control, TS_control = TS_control, 
       nsubsets = nsubsets, subset_rule = subset_rule,
       soften = soften, quiet = quiet)
}


