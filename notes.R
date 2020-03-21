# work on making formulas be able to be repeatedly prepped
# its not really right to use gamma as a name for ilr 
# lots of tidying to do!

devtools::load_all()

   data(rodents)

data <- rodents
topics = 2 
reps = 2
formulas = ~ 1
nchangepoints = 0:2
timename = "newmoon"
weights = TRUE
control = list()
control <- do.call("LDA_TS_control", control)
data <- conform_data(data = data, control = control)
LDAs <- LDA(data = data, topics = topics, reps = reps, 
              control = control$LDA_control)
control$TS_control$method_args$nit <- 100
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
