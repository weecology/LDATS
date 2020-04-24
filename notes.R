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
allow TS to alternatively have a data input and then run through identityLDA
turn the vcov code in est_regressors into functions/methpds  

something weird with ren's data having a -Inf lm logLik
ugh looks like the vcv can sometimes be singular
but really a broader issue is that most of the cp locations give -Infs
 which mean it's not a good fitting model, so how to address?
  at least make sure all initial cp locations are not -Inf, otherwise breaks!
hmmm yeah, think on this---how best to set up a way to handle models that
are just not going to fit
also we should probably allow the inclusion of starting values for params,
esp change points here....that's the issue

have a stop-gap patch that looks like it's working, but finding the root 
cause of that issue will be helpful in the long run

devtools::load_all()
   data(rodents)
LDAs <- LDA(data = rodents, topics = 2, replicates = c(1))

rm(list=ls())
LDAs <- LDA(data = rodents, topics = 4, replicates = c(1))

TSs <- TS(LDAs = LDAs, formulas = ~ 1, nchangepoints = 1, 
          timename = "newmoon", weights = TRUE,
          control = list(response = simplex_TS, 
                         response_args = list(control = 
                                              list(transformation = "ilr")),
                         method_args = list(control = list(nit = 100))))



TSs <- TS(LDAs = LDAs, formulas = ~ 1, nchangepoints = 0, 
          timename = "newmoon", #weights = TRUE,
          control = list(response = simplex_TS, 
                         response_args = list(control = 
                                              list(transformation = "clr")),
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