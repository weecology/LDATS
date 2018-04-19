data(rodents)
lda_data <- dplyr::select(rodents, -c(newmoon, date, plots, traps))
ts_data <- data.frame("time" = data.frame(rodents)[ , "newmoon"])

r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 4)

r_LDATS <- LDATS::LDA_TS(lda_data, ts_data, formula = c("1"),
                         ntopics = 2:5, nseeds = 2, ncores = 4, nit = 100)

document_term_matrix = lda_data
document_covariate_matrix = ts_data
formula = "1"
 nchangepoints = 1

MTS_set
data <- mtss
formula 
nchangepoints
weights


oo <- LDATS::MTS_set(data = mtss, formula = c("1", "time"), nchangepoints = 1, 
           weights = weights, nit = 10)


formula <- "1"

nchangepoints <- 1
nit <- 10
magnitude = 12
data <- data[[1]]

LDATS::prep_changepts(data, formula, ntemps, nchangepoints, weights)



