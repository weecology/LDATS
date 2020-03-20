# work on making formulas be able to be repeatedly prepped

# in the midst of the TS classic ldats...looks like starting to work, need
# to double back through tho

# get memo into a better location

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

control <- list()
   lda <- LDA(rodents, 3, 1, list(quiet = TRUE))

lda

   dct <- rodents$document_covariate_table
   dct$gamma <- lda[[1]]@gamma
   weights <- document_weights(dtt)
   check_multinom_TS_inputs(dct, timename = "newmoon")
   mts <- multinom_TS(dct, formula = gamma ~ 1, changepoints = c(20,50),
                      timename = "newmoon", weights = weights) 
 formula = gamma ~ 1
  data <- dct
  fit <- multinom(formula, data, weights, subset = in_chunk, trace = FALSE,
                  decay = control$lambda)

  basis <- "ilr"
  props <- data[ , grepl("gamma", colnames(data))]
  data$coords <- coda.base::coordinates(props, basis)
  formula <- as.formula(gsub("gamma", "coords", deparse(formula)))
  ff <-  lm(formula, data, subset = in_chunk, weights)
  fv_nrows <- sum(in_chunk)
  fv_ncols <- sum(grepl("gamma", colnames(data)))
  ff_fv <- as.matrix(ff$fitted.values, nrows = fv_nrows, ncols = fv_ncols)
  comps <- coda.base::composition(ff$fitted.values, basis)
