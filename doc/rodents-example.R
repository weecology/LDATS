## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE------------------------------------------------------
library(LDATS)
vers <- packageVersion("LDATS")
today <- Sys.Date()

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("weecology/LDATS")
#  library(LDATS)

## ------------------------------------------------------------------------
data(rodents)
head(rodents$document_term_table, 10)
head(rodents$document_covariate_table, 10)

## ----lda_set, eval =F----------------------------------------------------
#  lda_model_set <- LDA_set(document_term_table = rodents$document_term_table,
#                           topics = c(2:5),
#                           nseeds = 10,
#                           control = list(quiet = TRUE))
#  

## ----lda set not quiet, eval =F------------------------------------------
#  lda_model_set2 <- LDA_set(document_term_table = rodents$document_term_table,
#                           topics = c(2:3),
#                           nseeds = 2)

## ----load lda model set, include = F-------------------------------------
load(here::here('vignettes', 'rodents-example-files', 'lda_model_set.Rds'))
rm(lda_model_set2)

## ----select LDA----------------------------------------------------------
selected_lda_model <- select_LDA(lda_model_set)

## ----LDA results---------------------------------------------------------
# Number of topics:

selected_lda_model[[1]]@k

# Topic composition of communities at each time step
# Columns are topics; rows are time steps
head(selected_lda_model[[1]]@gamma)


## ----plot lda, fig.width=7, fig.height=6---------------------------------
plot(selected_lda_model[[1]])


## ----ts on lda, eval = F-------------------------------------------------
#  changepoint_models <- TS_on_LDA(LDA_models = selected_lda_model,
#                                  document_covariate_table = rodents$document_covariate_table,
#                                  formulas = ~ sin_year + cos_year,
#                                  nchangepoints = c(0:1),
#                                  timename = "newmoon",
#                                  weights = document_weights(rodents$document_term_table),
#                                  control = list(nit = 1000))
#  

## ----reload ts, include = F----------------------------------------------
load(here::here('vignettes', 'rodents-example-files', 'changepoint_models.Rds'))

## ----select ts-----------------------------------------------------------
selected_changepoint_model <- select_TS(changepoint_models)

## ----cpt results---------------------------------------------------------
# Number of changepoints
selected_changepoint_model$nchangepoints

# Summary of timesteps (newmoon values) for each changepoint
selected_changepoint_model$rho_summary

# Raw estimates for timesteps for each changepoint
# Changepoints are columns
head(selected_changepoint_model$rhos)


## ----plot cpt, fig.width=7, fig.height=6---------------------------------
plot(selected_changepoint_model)

## ----lda_ts, eval = F----------------------------------------------------
#  lda_ts_results <- LDA_TS(data = rodents,
#                           nseeds = 10,
#                           topics = 2:5,
#                           formulas = ~ sin_year + cos_year,
#                           nchangepoints= 0:1,
#                           timename = "newmoon",
#                           control = list(nit = 1000))

## ----load ldats results, include = F-------------------------------------
load(here::here('vignettes', 'rodents-example-files', 'lda_ts_results.Rds'))

## ----LDA_TS results------------------------------------------------------
names(lda_ts_results)

# Number of topics
lda_ts_results$`Selected LDA model`$k@k

# Number of changepoints
lda_ts_results$`Selected TS model`$nchangepoints

# Summary of changepoint locations
lda_ts_results$`Selected TS model`$rho_summary

## ----plot LDA_TS results, fig.height = 16, fig.width = 7, echo = F-------
plot(lda_ts_results)

