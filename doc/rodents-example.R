## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE------------------------------------------------------
library(LDATS)
vers <- packageVersion("LDATS")
today <- Sys.Date()

## ----download files, include = FALSE-------------------------------------
  vignette_files <- tempdir()
  dir.create(file.path(vignette_files, "output"), showWarnings = FALSE)
  github_path <- "https://github.com/weecology/LDATS-replications/raw/master/output/"
  files_to_download <- c("rodents_example_lda_model_set.RDS", "rodents_example_ts_model_set.RDS", 
                         "rodents_example_lda_ts_model_set.RDS")
  
  for (file in files_to_download)  {
    download.file(url = paste0(github_path, file),
                  destfile = file.path(vignette_files, "output", file), 
                  mode = "wb")
  }

## ---- eval = FALSE-------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("weecology/LDATS")

## ------------------------------------------------------------------------
data(rodents)
head(rodents$document_term_table, 10)
head(rodents$document_covariate_table, 10)

## ----lda_set, eval = FALSE-----------------------------------------------
#  lda_model_set <- LDA(data = rodents, topics = 2:5, replicates = 10,
#                       control = list(quiet = TRUE))

## ----save lda model set, include = FALSE, eval = FALSE-------------------
#  saveRDS(lda_model_set, file.path(vignette_files, "output", "rodents_example_lda_model_set.RDS"))

## ----lda set not quiet, eval = FALSE-------------------------------------
#  lda_model_set2 <- LDA(data = rodents, topics = c(2:3), replicates = 2)

## ----load lda model set, include = FALSE---------------------------------
lda_model_set <- readRDS(file.path(vignette_files, "output", "rodents_example_lda_model_set.Rds"))
rm(lda_model_set2)

## ----select LDA----------------------------------------------------------
selected_lda_model <- select_LDA(lda_model_set$LDAs)[[1]]

## ----LDA results---------------------------------------------------------
# Number of topics:

selected_lda_model$topics

# Topic composition of communities at each time step
# Columns are topics; rows are time steps
head(selected_lda_model$document_topic_table)

## ----plot lda, fig.width = 7, fig.height = 6-----------------------------
plot(selected_lda_model)

## ----ts set, eval = FALSE------------------------------------------------
#  ts_model_set <- TS(LDAs = lda_model_set,
#                     formulas = ~ sin_year + cos_year,
#                     nchangepoints = 0:1,
#                     timename = "newmoon",
#                     weights = TRUE,
#                     control = list(method_args =
#                             list(control = ldats_classic_control(nit = 1000))))

## ----save ts model set, include = FALSE, eval = FALSE--------------------
#  saveRDS(ts_model_set, file.path(vignette_files, "output", "rodents_example_ts_model_set.RDS"))

## ----load ts model set, include = FALSE----------------------------------
ts_model_set <- readRDS(file.path(vignette_files, "output", "rodents_example_ts_model_set.RDS"))

## ----select ts-----------------------------------------------------------
selected_changepoint_model <- select_TS(ts_model_set$TSs)[[1]]

## ----cpt results---------------------------------------------------------
# Number of changepoints
selected_changepoint_model$nchangepoints

# Summary of timesteps (newmoon values) for each changepoint
selected_changepoint_model$rho_summary

# Raw estimates for timesteps for each changepoint
# Changepoints are columns
head(selected_changepoint_model$focal_rhos)


## ----plot cpt, fig.width=7, fig.height=6---------------------------------
plot(selected_changepoint_model)

## ----lda_ts, eval = FALSE------------------------------------------------
#  lda_ts_results <- LDA_TS(data = rodents,
#                           replicates = 10,
#                           topics = 2:5,
#                           formulas = ~ sin_year + cos_year,
#                           nchangepoints= 0:1,
#                           timename = "newmoon",
#                           control = list(TS_method_args =
#                             list(control = ldats_classic_control(nit = 1000))))

## ----save lda ts model set, include = FALSE, eval = FALSE----------------
#  saveRDS(lda_ts_results, file.path(vignette_files, "output", "rodents_example_lda_ts_model_set.RDS"))

## ----load lda ts model set, include = FALSE------------------------------
lda_ts_results <- readRDS(file.path(vignette_files, "output", "rodents_example_lda_ts_model_set.RDS"))

## ----LDA_TS results------------------------------------------------------
names(lda_ts_results)

# Number of topics
lda_ts_results$"LDA models"$selected_LDAs[[1]]$topics

# Number of changepoints
lda_ts_results$"TS models"$selected_TSs[[1]]$nchangepoints

# Summary of changepoint locations
lda_ts_results$"TS models"$selected_TSs[[1]]$rho_summary

## ----plot LDA_TS results, fig.height = 16, fig.width = 7, echo = F-------
plot(lda_ts_results)

