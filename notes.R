# work on making formulas be able to be repeatedly prepped
# its not really right to use gamma as a name for ilr 
# lots of tidying to do!


known to dos
pipeline generalization and doc streamline
check functions
utilities
plotting functions
simulate
generalize ilr to simplical to allow for other transforms


devtools::load_all()

   data(rodents)

data <- rodents
topics = 2 
replicates = 2
formulas = ~ 1
nchangepoints = 1
timename = "newmoon"
weights = TRUE
control = list()
LDAs <- LDA(data = data, topics = topics, replicates = replicates, 
              control = list())
control <- LDA_TS_control()
control$TS_control$method_args$control$nit <- 100
TSs <- TS(LDAs = LDAs, data = data, formulas = formulas, 
          nchangepoints = nchangepoints, timename = timename, 
          weights = weights, control = control$TS_control)


now LDA could be any of the linquistic decomposition analyses or whatever
including any of a number of "LDA" functions or models, so i'm revoking
the importing of LDA from topicmodels but keeping topicmodels imported
to allow calling of topicmodels::LDA from inside LDA (the default!)


seed is rep

introduction of soften logical variable
designed to help soften errors in pipelines through wrapping in tryCatch

LDA_set is now LDA and TS_on_LDA is now TS

standardized output from LDA and TS functions


# LDA is now functionally what LDA_set was before
# "reps" replaces "nseeds" but otherwise things are the same generally at api
# the control list is expanded to now include a base LDA modeling function
# and then arguments for all three functions, it also gets the subset info
# because it may be needed, and in addition to quiet a boolean "soften"
# which softens errors during model running
# prep_LDA_models combines the data (conforming if needed), topics, 
# reps, and subsets to prep the list of model inputs, that will actually be
# appended to with the results, which are done at the top with LDA_call.
# LDA_call prints the message, and then preps the input arguments as used
# in a do.call on the function given in the control list. there are
# examples of both a generalized topicmodels LDA function wrapper (works
# with Gibbs sampling approach too!) and an "identity" 1-topic function.
# these functions allow for pre- and post-processing around the main model
# function without having to impact the general LDA_call
# select_LDA does a general measuring and then selecting, which now allows
# for arguments to be passed in!


# data can come into LDA_TS LDA TS in a variety of forms, and depending on 
# usages, might take a variety of different forms
# the purpose of this function is to generalize and extract the code used
# to shuddle between data formats from functions / replace with a single line
# it's still a work in progress and needs more extensive usage exploration,
# as it's going to be a workhorse function. 
currently its not saving the test/train split explicitly,
# just implicitly via the data encoding that exists. we should probably
# shore this up a bit more for sure.
# also this function is big and modularized a good degree already...it could
# get chunked into subfunctions
# there are functions for basic leave p out cross validation, including
# both systematic and random approaches


LDA_TS

  LDA_TS_control
    LDA_control
    TS_control

  LDA
    prepare_LDA
      LDA_control
      conform_data
    run_LDA
      LDA_call (replicated for each model)
        LDA_msg
        LDA$control$model with LDA$control$model_args
          topicmodels::LDA
    package_LDA
      select_LDA
        LDA$control$selector with LDA$control$selector_args
        measure_LDA (replicated for each model)
          LDA$control$measurer with LDA$control$measurer_args

  TS
    prepare_TS
      TS_control
    run_TS
      TS_call (replicated for each model)
        TS_msg
        TS$control$model with TS$control$model_args
          sequential_TS_control
          est_changepoints
          est_regressors
          package_sequential_TS
    package_TS
      select_TS
        TS$control$selector with TS$control$selector_args
        measure_TS (replicated for each model)
          TS$control$measurer with TS$control$measurer_args



