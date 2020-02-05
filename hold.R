   data(rodents)
   dtt <- rodents$document_term_table
   lda <- LDA_set(dtt, 3, 1, list(quiet = TRUE))
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
