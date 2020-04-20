to-do

predict
tidying check functions and general usage
examples
tests
tempering
vignettes
  

devtools::load_all()
   data(rodents)
LDAs <- LDA(data = rodents, topics = 2:3, replicates = c(1,4))


TS(LDAs = LDAs, formulas = ~ newmoon, nchangepoints = 1, 
          timename = "newmoon", weights = TRUE,
          control = list(response = simplex_TS, 
                         method_args = 
                         list(control = ldats_classic_control(nit = 100)),
                                response_args = 
                           list(control = 
                            simplex_TS_control(transformation = "ilr"))))

TS(LDAs = LDAs, formulas = ~ newmoon, nchangepoints = 1, 
          timename = "newmoon", weights = TRUE,
          control = list(method_args = 
                         list(control = ldats_classic_control(nit = 100))))



plot(LDAs)
plot(TSs)
plot(TSs, plot_type="diagnostic")
plot(LDATSs)


not sure if needed:

time_order_data <- function(x, timename = "time"){
  time_order <- order(x[ , timename])
  x[time_order , ]
}