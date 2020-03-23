predict
tidying check functions and general usage
examples
tests
what is up with lda_ts_controls...
  ...oh is that how its always been because of the splitting?
tempering
vignettes

devtools::load_all()
   data(rodents)
LDAs <- LDA(data = rodents, topics = 1:2, replicates = c(1, 4))
TSs <- TS(LDAs = LDAs, formulas = ~ 1, nchangepoints = 1, 
          timename = "newmoon", weights = TRUE,
          control = list(method_args = 
                         list(control = ldats_classic_control(nit = 100))))
LDATSs <- LDA_TS(data = rodents, topics = 1:2, replicates = c(1, 4),
                 formulas = ~ 1, nchangepoints = 1, timename = "newmoon",
                 weights = TRUE,
                 control = list(TS_method_args = 
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