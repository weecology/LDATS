plotting functions
  -generalize pred_gamma_TS_plot
simulate
  -just generalize TS 
predict
tidying check functions
tempering
examples
tests

what is up with lda_ts_controls...
  ...oh thats how its always been because of the splitting?
where do the ts_fit functions go? ldats_classic?


devtools::load_all()
   data(rodents)
data <- rodents
topics = 1:2 
replicates = c(2,5)
formulas = ~ 1
nchangepoints = 1
timename = "newmoon"
weights = TRUE
control = 
LDAs <- LDA(data = data, topics = topics, replicates = replicates, 
              control = list())
control$TS_control$method_args$control$nit <- 100
TSs <- TS(LDAs = LDAs, formulas = formulas, 
          nchangepoints = nchangepoints, timename = timename, 
          weights = weights, control = 
list(method_args = list(control = ldats_classic_control(nit = 100))))


LDATSs <- LDA_TS(data = data, topics = topics, replicates = replicates,
                 formulas = formulas, 
          nchangepoints = nchangepoints, timename = timename, 
          weights = weights, control = list(TS_method_args = 
                          list(control = ldats_classic_control(nit=100))))

names(LDATSs)

cc<-TS_control(method_args = list(control = ldats_classic_control(nit = 100)))
$method_args$control$nit

cc<-LDA_TS_control(TS_method_args = list(control = ldats_classic_control(nit = 100)))

names(cc)
names(formals(TS_control))

$TS_control$method_args$control$nit



not sure if needed:

time_order_data <- function(x, timename = "time"){
  time_order <- order(x[ , timename])
  x[time_order , ]
}