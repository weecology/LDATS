plotting functions
  -TS and LDA_TS
simulate
  -just generalize TS 
predict
tidying check functions
tempering
examples

devtools::load_all()
   data(rodents)
data <- rodents
topics = 1:2 
replicates = c(2,5)
formulas = ~ 1
nchangepoints = 1
timename = "newmoon"
weights = TRUE
control = LDA_TS_control()
LDAs <- LDA(data = data, topics = topics, replicates = replicates, 
              control = control$LDA_control)
control$TS_control$method_args$control$nit <- 100
TSs <- TS(LDAs = LDAs, formulas = formulas, 
          nchangepoints = nchangepoints, timename = timename, 
          weights = weights, control = control$TS_control)

not sure if needed:

time_order_data <- function(x, timename = "time"){
  time_order <- order(x[ , timename])
  x[time_order , ]
}