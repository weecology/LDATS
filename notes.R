to-do

predict
tidying check functions and general usage
examples
tests
softcall message handling
tempering
vignettes
   data 
   simplex response
  

devtools::load_all()
   data(rodents)
LDAs <- LDA(data = rodents, topics = 2:3, replicates = c(1,4))


TSs <- TS(LDAs = LDAs, formulas = ~ newmoon, nchangepoints = 1, 
          timename = "newmoon", weights = TRUE,
          control = list(response = simplex_TS, 
                         response_args = list(control = 
                                              list(transformation = "ilr")),
                         method_args = list(control = list(nit = 100))))

TSs <- TS(LDAs = LDAs, formulas = ~ newmoon, nchangepoints = 1, 
          timename = "newmoon", weights = TRUE,
          control = list(method_args = list(control = list(nit = 100))))

plot(LDAs)
plot(TSs)
plot(TSs, plot_type="diagnostic")
plot(LDATSs)


not sure if needed:

time_order_data <- function(x, timename = "time"){
  time_order <- order(x[ , timename])
  x[time_order , ]
}