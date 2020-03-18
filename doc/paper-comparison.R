params <-
list(run_models = FALSE)

## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE-----------------------------------------------------------
vers <- packageVersion("LDATS")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("weecology/LDATS")

## -----------------------------------------------------------------------------
library(LDATS)
set.seed(42)
nseeds <- 200
nit <- 10000

## ---- eval = params$run_models------------------------------------------------
#  install.packages(c("dplyr", "gridExtra", "multipanel", "RColorBrewer", "RCurl", "reshape2"))

## ---- eval = FALSE------------------------------------------------------------
#  rmarkdown::render("paper-comparison.Rmd", params = list(run_models = TRUE))

## ----set download location----------------------------------------------------
vignette_files <- tempdir()

## ----download scripts and data------------------------------------------------
test_file <- file.path(vignette_files, "scripts", "rodent_LDA_analysis.r")

if (!file.exists(test_file)){

  # scripts
  dir.create(file.path(vignette_files, "scripts"), showWarnings = FALSE)
  github_path <- "https://raw.githubusercontent.com/weecology/LDATS-replications/master/scripts/"
  files_to_download <- c("rodent_LDA_analysis.r", "rodent_data_for_LDA.r", 
                         "AIC_model_selection.R", "changepointmodel.r", 
                         "LDA-distance.R", "LDA_figure_scripts.R")
  
  for (file in files_to_download)  {
    download.file(url = paste0(github_path, file),
                  destfile = file.path(vignette_files, "scripts", file))
  }

    
  # data
  dir.create(file.path(vignette_files, "data"), showWarnings = FALSE)
  github_path <- "https://raw.githubusercontent.com/weecology/LDATS-replications/master/data/"
  files_to_download <- c("moon_dates.csv", "Portal_rodent_trapping.csv", 
                         "Rodent_table_dat.csv", "paper_dat.csv",
                         "paper_dates.csv", "paper_covariates.csv")
  
  for (file in files_to_download)  {
    download.file(url = paste0(github_path, file),
                  destfile = file.path(vignette_files, "data", file))
  }
}

## ----download pre-generated model outputs-------------------------------------
test_file <- file.path(vignette_files, "output", "ldats_ldamodel.RDS")

if (!file.exists(test_file)){

  dir.create(file.path(vignette_files, "output"), showWarnings = FALSE)
  github_path <- "https://raw.githubusercontent.com/weecology/LDATS-replications/master/output/"
  files_to_download <- c("ldats_ldamodel.RDS", "paper_ldamodel.RDS", 
                         "ldats_ldats.RDS", "ldats_paper.RDS", 
                         "paper_ldats.RDS", "paper_paper.RDS", 
                         "ldats_rodents_adjusted.RDS", "rodents.RDS",
                         "ldats_paper_cpt.RDS", "ldats_paper_cpt_dates.RDS",
                         "ldats_ldats_cpt.RDS", "ldats_ldats_cpt_dates.RDS",
                         "paper_paper_cpt.RDS", "paper_paper_cpt_dates.RDS",
                         "paper_ldats_cpt.RDS", "paper_ldats_cpt_dates.RDS",
                         "annual_hist.RDS", "cpt_dates.RDS",
                         "lda_distances.png", "paper_paper_cpt_plot.png",
                         "ldats_paper_cpt_plot.png")

  for (file in files_to_download){
    download.file(url = paste0(github_path, file),
                  destfile = file.path(vignette_files, "output", file), 
                  mode = "wb")
  }
}

## ----LDATS data---------------------------------------------------------------
data(rodents)

head(rodents[[1]])

## ----Paper data, eval = params$run_models-------------------------------------
#  # parameters for subsetting the full Portal rodents data
#  periods <- 1:436
#  control_plots <- c(2, 4, 8, 11, 12, 14, 17, 22)
#  species_list <- c("BA", "DM", "DO", "DS", "NA", "OL", "OT", "PB", "PE", "PF",
#                    "PH", "PI", "PL", "PM", "PP", "RF", "RM", "RO", "SF", "SH", "SO")
#  
#  source(file.path(vignette_files, "scripts", "rodent_data_for_LDA.r"))
#  
#  # assemble `paper_dat`, the data from Christensen et al. 2018
#  paper_dat <- create_rodent_table(period_first = min(periods),
#                                   period_last = max(periods),
#                                   selected_plots = control_plots,
#                                   selected_species = species_list)
#  
#  # assemble `paper_covariates`, the associated dates and covariate data
#  moondat <- read.csv(file.path(vignette_files, "data", "moon_dates.csv"), stringsAsFactors = F)
#  
#  paper_dates <- moondat %>%
#    dplyr::filter(period %>% dplyr::between(min(periods), max(periods))) %>%
#    dplyr::pull(censusdate) %>%
#    as.Date()
#  
#  paper_covariates <- data.frame(
#    index = seq_along(paper_dates),
#    date = paper_dates,
#    year_continuous = lubridate::decimal_date(paper_dates)) %>%
#    dplyr::mutate(
#      sin_year = sin(year_continuous * 2 * pi),
#      cos_year = cos(year_continuous * 2 * pi)
#    )

## ----Paper data2, eval = !params$run_models, include = FALSE------------------
  moondat <- read.csv(file.path(vignette_files, "data", "moon_dates.csv"), stringsAsFactors = FALSE)
  paper_dat <- read.csv(file.path(vignette_files, "data", "paper_dat.csv"), stringsAsFactors = FALSE)
  paper_dates <- read.csv(file.path(vignette_files, "data", "paper_dates.csv"), stringsAsFactors = FALSE)
  paper_covariates <- read.csv(file.path(vignette_files, "data", "paper_covariates.csv"), stringsAsFactors = FALSE)

## ----rodent data comparison---------------------------------------------------
compare <- rodents[[1]] == paper_dat

length(which(rowSums(compare) < ncol(compare)))

## ----Data adjustment eval, eval = FALSE---------------------------------------
#    # retrieve data on number of plots trapped per month
#    trap_table = read.csv('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv')
#    trap_table_controls = filter(trap_table, plot %in% selected_plots)
#    nplots_controls = aggregate(trap_table_controls$sampled,by=list(period = trap_table_controls$period),FUN=sum)
#  
#    # adjust species counts by number of plots trapped that month
#    r_table_adjusted = as.data.frame.matrix(r_table)
#    for (n in 1:436) {
#      #divide by number of control plots actually trapped (should be 8) and multiply by 8 to estimate captures as if all plots were trapped
#      r_table_adjusted[n,] = round(r_table_adjusted[n,]/nplots_controls$x[n]*8)
#    }

## ----adjust LDATS data after Christensen et al, eval = params$run_models------
#  # get the trapping effort for each sample
#  trap_table <- read.csv(file.path(vignette_files, "data", "Portal_rodent_trapping.csv"))
#  trap_table_controls <- dplyr::filter(trap_table, plot %in% control_plots)
#  nplots_controls <- aggregate(trap_table_controls$sampled,
#                              by = list(period = trap_table_controls$period),
#                              FUN = sum)
#  
#  # adjust species counts by number of plots trapped that month
#  #   divide by number of control plots actually trapped (should be 8) and
#  #   multiply by 8 to estimate captures as if all plots were trapped
#  ldats_rodents_adjusted <- as.data.frame.matrix(rodents[[1]])
#  ldats_rodents_adjusted[periods, ] <- round(ldats_rodents_adjusted[periods, ] / nplots_controls$x[periods] * 8)

## ----eval = params$run_models, include = FALSE--------------------------------
#  saveRDS(ldats_rodents_adjusted, file = file.path(vignette_files, "output", "ldats_rodents_adjusted.RDS"))

## ----eval = !params$run_models, include = FALSE-------------------------------
ldats_rodents_adjusted <- readRDS(file.path(vignette_files, "output", "ldats_rodents_adjusted.RDS"))

## ----dataset comparisons, eval = params$run_models----------------------------
#  compare_raw <- rodents[[1]] == ldats_rodents_adjusted
#  length(which(rowSums(compare_raw) < ncol(compare_raw)))
#  
#  compare_adjusted <- ldats_rodents_adjusted == paper_dat
#  length(which(rowSums(compare_adjusted) < ncol(compare_adjusted)))

## ----switch to adjusted rodents-----------------------------------------------
rodents[[1]] <- paper_dat

## ----show covariate table-----------------------------------------------------
head(rodents$document_covariate_table)

## ----eval = params$run_models, include = FALSE--------------------------------
#  "%>%" <- dplyr::"%>%"

## ----add dates to covariate table, eval = params$run_models-------------------
#  new_cov_table <- dplyr::left_join(rodents$document_covariate_table,
#                                    dplyr::select(moondat, newmoonnumber, censusdate),
#                                    by = c("newmoon" = "newmoonnumber")) %>%
#                                    dplyr::rename(date = censusdate)
#  
#  rodents$document_covariate_table <- new_cov_table

## ----eval = params$run_models, include = FALSE--------------------------------
#  saveRDS(rodents, file = file.path(vignette_files, "output", "rodents.RDS"))

## ----LDATS LDAs, eval = params$run_models-------------------------------------
#  ldats_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table,
#                               topics = 2:6, nseeds = nseeds)
#  ldats_ldamodel <- LDATS::select_LDA(LDA_models = ldats_ldas)[[1]]
#  
#  saveRDS(ldats_ldamodel, file = file.path(vignette_files, "ldats_ldamodel.RDS"))

## ----paper LDAs, eval = params$run_models-------------------------------------
#  source(file.path(vignette_files, "scripts", "AIC_model_selection.R"))
#  source(file.path(vignette_files, "scripts", "LDA-distance.R"))
#  
#  # Some of the functions require the data to be stored in the `dat` variable
#  dat <- paper_dat
#  
#  # Fit a bunch of LDA models with different seeds
#  # Only use even numbers for seeds because consecutive seeds give identical results
#  seeds <- 2 * seq(nseeds)
#  
#  # repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
#  best_ntopic <- repeat_VEM(paper_dat,
#                            seeds,
#                            topic_min = 2,
#                            topic_max = 6)
#  hist(best_ntopic$k, breaks = seq(from = 0.5, to = 9.5),
#       xlab = "best # of topics", main = "")
#  
#  # 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
#  # ==================================================================
#  # get the best 100 seeds where 4 topics was the best LDA model
#  seeds_4topics <- best_ntopic %>%
#    filter(k == 4) %>%
#    arrange(aic) %>%
#    head(min(100, nseeds)) %>%
#    pull(SEED)
#  
#  # choose seed with highest log likelihood for all following analyses
#  #    (also produces plot of community composition for "best" run compared to "worst")
#  
#  png(file.path(vignette_files, "output", "lda_distances.png"), width = 800, height = 400)
#  dat <- paper_dat # calculate_LDA_distance has some required named variables
#  best_seed <- calculate_LDA_distance(paper_dat, seeds_4topics)
#  dev.off()
#  mean_dist <- unlist(best_seed)[2]
#  max_dist <- unlist(best_seed)[3]
#  
#  # ==================================================================
#  # 3. run LDA model
#  # ==================================================================
#  ntopics <- 4
#  SEED <- unlist(best_seed)[1]  # For the paper, use seed 206
#  ldamodel <- LDA(paper_dat, ntopics, control = list(seed = SEED), method = "VEM")
#  
#  saveRDS(ldamodel, file = file.path(vignette_files, "paper_ldamodel.RDS"))

## -----------------------------------------------------------------------------
knitr::include_graphics(file.path(vignette_files, "output", "lda_distances.png"))

## -----------------------------------------------------------------------------
ldamodel <- readRDS(file.path(vignette_files, "output", "paper_ldamodel.RDS"))
ldats_ldamodel <- readRDS(file.path(vignette_files, "output", "ldats_ldamodel.RDS"))

## ----plot paper LDA, fig.width = 7, fig.height = 6----------------------------
plot(ldamodel, cols = NULL, option = "D")

## ----plot LDATS LDA, fig.width = 7, fig.height = 6----------------------------
plot(ldats_ldamodel, cols = NULL, option = "D")

## ----paper changepoint models, eval = params$run_models-----------------------
#  #### Run changepoint ####
#  source(file.path(vignette_files, "scripts", "changepointmodel.r"))
#  
#  find_changepoints <- function(lda_model, paper_covariates, n_changepoints = 1:6){
#    # set up parameters for model
#    x <- dplyr::select(paper_covariates,
#                       year_continuous,
#                       sin_year,
#                       cos_year)
#  
#    # run models with 1, 2, 3, 4, 5, 6 changepoints
#    cpt_results <- data.frame(n_changepoints = n_changepoints)
#    cpt_results$cpt_model <- lapply(cpt_results$n_changepoints,
#                                    function(n_changepoints){
#                                      changepoint_model(lda_model, x, n_changepoints, maxit = nit,
#                                                        weights = rep(1, NROW(x)))
#                                    })
#    return(cpt_results)
#  }
#  
#  # Among a selection of models with different # of changepoints,
#  #   - compute AIC
#  #   - select the model with the best AIC
#  #   - get the posterior distributions for the changepoints
#  select_cpt_model <- function(cpt_results, ntopics){
#    # compute log likelihood as the mean deviance
#    cpt_results$mean_deviances <- vapply(cpt_results$cpt_model,
#                                         function(cpt_model) {mean(cpt_model$saved_lls)},
#                                         0)
#  
#    # compute AIC = ( -2 * log likelihood) + 2 * (#parameters)
#    cpt_results$AIC <- cpt_results$mean_deviances * -2 +
#      2 * (3 * (ntopics - 1) * (cpt_results$n_changepoints + 1) +
#             (cpt_results$n_changepoints))
#  
#    # select the best model
#    cpt <- cpt_results$cpt_model[[which.min(cpt_results$AIC)]]
#    return(cpt)
#  }
#  
#  # transform the output from `compute_cpt` and match up the time indices with
#  #   dates from the original data
#  get_dates <- function(cpt, covariates = paper_covariates){
#    cpt$saved[,1,] %>%
#      t() %>%
#      as.data.frame() %>%
#      reshape::melt() %>%
#      dplyr::left_join(covariates, by = c("value" = "index"))
#  }

## ----save annual_hist, include = FALSE, eval = params$run_models--------------
#  saveRDS(annual_hist, file = file.path(vignette_files, "output", "annual_hist.RDS"))

## ----run LDATS LDA and paper cpt, eval = params$run_models--------------------
#  ldats_paper_results <- find_changepoints(ldats_ldamodel, paper_covariates)
#  
#  saveRDS(ldats_paper_results, file = file.path(vignette_files, "output", "ldats_paper.RDS"))

## ----compute changepoints for LDATS LDA and paper cpt, eval = params$run_models----
#  ldats_paper_results <- readRDS(file.path(vignette_files, "output", "ldats_paper.RDS"))
#  
#  ldats_paper_cpt <- select_cpt_model(ldats_paper_results,
#                                      ntopics = ldats_ldamodel@k)
#  ldats_paper_cpt_dates <- get_dates(ldats_paper_cpt)

## ----include = FALSE, eval = params$run_models--------------------------------
#  saveRDS(ldats_paper_cpt, file = file.path(vignette_files, "output", "ldats_paper_cpt.RDS"))
#  saveRDS(ldats_paper_cpt_dates, file = file.path(vignette_files, "output", "ldats_paper_cpt_dates.RDS"))

## ----run paper LDA and paper cpt, eval = params$run_models--------------------
#  paper_paper_results <- find_changepoints(ldamodel, paper_covariates)
#  
#  saveRDS(paper_paper_results, file = file.path(vignette_files, "paper_paper.RDS"))

## ----compute changepoints for paper LDA and paper cpt, eval = params$run_models----
#  paper_paper_results <- readRDS(file.path(vignette_files, "output", "paper_paper.RDS"))
#  
#  paper_paper_cpt <- select_cpt_model(paper_paper_results,
#                                      ntopics = ldamodel@k)
#  paper_paper_cpt_dates <- get_dates(ldats_paper_cpt)

## ----include = FALSE, eval = params$run_models--------------------------------
#  saveRDS(paper_paper_cpt, file = file.path(vignette_files, "output", "paper_paper_cpt.RDS"))
#  saveRDS(paper_paper_cpt_dates, file = file.path(vignette_files, "output", "paper_paper_cpt_dates.RDS"))

## ----run LDATS LDA and LDATS cpt, eval = params$run_models--------------------
#  ldats_ldats_results <- TS_on_LDA(LDA_models = ldats_ldamodel,
#                                   document_covariate_table = rodents$document_covariate_table,
#                                   formulas = ~ sin_year + cos_year,
#                                   nchangepoints = 1:6,
#                                   timename = "newmoon",
#                                   weights = NULL,
#                                   control = list(nit = nit))
#  
#  saveRDS(ldats_ldats_results, file = file.path(vignette_files, "output", "ldats_ldats.RDS"))

## ----construct lookup table for LDATS output for changepoint times, eval = params$run_models----
#  # make the full sequence of possible newmoon values
#  full_index <- seq(min(rodents$document_covariate_table$newmoon),
#                    max(rodents$document_covariate_table$newmoon))
#  
#  # generate a lookup table with dates for the newmoons, using `approx` to
#  #   linearly interpolate the missing values
#  ldats_dates <- approx(rodents$document_covariate_table$newmoon,
#                       as.Date(rodents$document_covariate_table$date),
#                       full_index) %>%
#    as.data.frame() %>%
#    mutate(index = x,
#           date = as.Date(y, origin = "1970-01-01")) %>%
#    select(index, date)

## ----compute changepoints for LDATS LDA and LDATS cpt, eval = params$run_models----
#  ldats_ldats_results <- readRDS(file.path(vignette_files, "output", "ldats_ldats.RDS"))
#  
#  ldats_ldats_cpt <- select_TS(ldats_ldats_results)
#  
#  ldats_ldats_cpt_dates <- ldats_ldats_cpt$rhos %>%
#    as.data.frame() %>%
#    reshape::melt() %>%
#    dplyr::left_join(ldats_dates, by = c("value" = "index"))

## ----include = FALSE, eval = params$run_models--------------------------------
#  saveRDS(ldats_ldats_cpt, file = file.path(vignette_files,  "output", "ldats_ldats_cpt.RDS"))
#  saveRDS(ldats_ldats_cpt_dates, file = file.path(vignette_files,  "output", "ldats_ldats_cpt_dates.RDS"))

## ----run paper LDA and LDATS cpt, eval = params$run_models--------------------
#  paper_ldats_results <- TS_on_LDA(LDA_models = ldamodel,
#                               document_covariate_table = rodents$document_covariate_table,
#                               formulas = ~ sin_year + cos_year,
#                               nchangepoints = 1:6,
#  
#                               timename = "newmoon",
#                               weights = NULL,
#                               control = list(nit = nit))
#  
#  
#  saveRDS(paper_ldats_results, file = file.path(vignette_files, "output", "paper_ldats.RDS"))

## ----select paper lda + ldats cpt, eval = params$run_models-------------------
#  paper_ldats_results <- readRDS(file.path(vignette_files, "output", "paper_ldats.RDS"))
#  
#  paper_ldats_cpt <- select_TS(paper_ldats_results)
#  
#  paper_ldats_cpt_dates <- paper_ldats_cpt$rhos %>%
#    as.data.frame() %>%
#    reshape::melt() %>%
#    dplyr::left_join(ldats_dates, by = c("value" = "index"))

## ----include = FALSE, eval = params$run_models--------------------------------
#  saveRDS(paper_ldats_cpt, file = file.path(vignette_files,  "output", "paper_ldats_cpt.RDS"))
#  saveRDS(paper_ldats_cpt_dates, file = file.path(vignette_files,  "output", "paper_ldats_cpt_dates.RDS"))

## ----eval = !params$run_models, include = FALSE-------------------------------
ldats_paper_cpt <- readRDS(file.path(vignette_files, "output", "ldats_paper_cpt.RDS"))
paper_paper_cpt <- readRDS(file.path(vignette_files, "output", "paper_paper_cpt.RDS"))
paper_ldats_cpt <- readRDS(file.path(vignette_files, "output", "paper_ldats_cpt.RDS"))
ldats_ldats_cpt <- readRDS(file.path(vignette_files, "output", "ldats_ldats_cpt.RDS"))
ldats_paper_cpt_dates <- readRDS(file.path(vignette_files, "output", "ldats_paper_cpt_dates.RDS"))
paper_paper_cpt_dates <- readRDS(file.path(vignette_files, "output", "paper_paper_cpt_dates.RDS"))
paper_ldats_cpt_dates <- readRDS(file.path(vignette_files, "output", "paper_ldats_cpt_dates.RDS"))
ldats_ldats_cpt_dates <- readRDS(file.path(vignette_files, "output", "ldats_ldats_cpt_dates.RDS"))

## -----------------------------------------------------------------------------
nlevels(ldats_paper_cpt_dates$variable)
nlevels(paper_paper_cpt_dates$variable)
nlevels(ldats_ldats_cpt_dates$variable)
nlevels(paper_ldats_cpt_dates$variable)

## ----plot paper LDA and LDATS cpts, fig.width = 7, fig.height = 6-------------
plot(paper_ldats_cpt)

## ----plot ldats LDA and LDATS cpt, fig.width = 7, fig.height = 6--------------
plot(ldats_ldats_cpt)

## ---- eval = !params$run_models-----------------------------------------------
annual_hist <- readRDS(file.path(vignette_files, "output", "annual_hist.RDS"))

## ----plot paper LDA and paper cpt, eval = params$run_models-------------------
#  paper_cpts <- find_changepoint_location(paper_paper_cpt)
#  ntopics <- ldamodel@k
#  
#  png(file.path(vignette_files, "output", "paper_paper_cpt_plot.png"), width = 800, height = 600)
#  get_ll_non_memoized_plot(ldamodel, paper_covariates, paper_cpts, make_plot = TRUE,
#                                             weights = rep(1, NROW(paper_covariates)))
#  dev.off()

## ----plot paper LDA and LDATS cpts2, fig.width = 7, fig.height = 6------------
paper_paper_hist <- annual_hist(paper_paper_cpt, paper_covariates$year_continuous)

## -----------------------------------------------------------------------------
knitr::include_graphics(file.path(vignette_files, "output", "paper_paper_cpt_plot.png"))

## ----plot LDATS lda and paper cpt, eval = params$run_models-------------------
#  ldats_cpts <- find_changepoint_location(ldats_paper_cpt)
#  ntopics <- ldats_ldamodel@k
#  
#  png(file.path(vignette_files, "output", "ldats_paper_cpt_plot.png"), width = 800, height = 600)
#  get_ll_non_memoized_plot(ldats_ldamodel, paper_covariates, ldats_cpts, make_plot = TRUE,
#                                             weights = rep(1, NROW(paper_covariates)))
#  dev.off()

## ----plot LDATS lda and paper cpt2, fig.width = 7, fig.height = 6-------------
ldats_paper_hist <- annual_hist(ldats_paper_cpt, paper_covariates$year_continuous)

## -----------------------------------------------------------------------------
knitr::include_graphics(file.path(vignette_files, "output", "ldats_paper_cpt_plot.png"))

## ----report cpt dates, include = FALSE, eval = params$run_models--------------
#  paper_paper_cpt_dates$date <- as.Date(paper_paper_cpt_dates$date)
#  ldats_paper_cpt_dates$date <- as.Date(ldats_paper_cpt_dates$date)
#  
#  cpt_dates <- dplyr::bind_rows("paperLDA_paperCPT" = paper_paper_cpt_dates,
#                                "ldatsLDA_paperCPT" = ldats_paper_cpt_dates,
#                                "ldatsLDA_ldatsCPT" = ldats_ldats_cpt_dates,
#                                "paperLDA_ldatsCPT" = paper_ldats_cpt_dates,
#                                .id = "analysis") %>%
#    dplyr::group_by(analysis, variable) %>%
#    dplyr::summarize(date = mean(date)) %>%
#    dplyr::ungroup() %>%
#    dplyr::rename(changepoint = variable) %>%
#    tidyr::spread(analysis, date)
#  
#  saveRDS(cpt_dates, file = file.path(vignette_files,  "output", "cpt_dates.RDS"))

## ----eval = !params$run_models, include = FALSE-------------------------------
cpt_dates <- readRDS(file.path(vignette_files, "output", "cpt_dates.RDS"))

## ----print cpt dates----------------------------------------------------------
knitr::kable(cpt_dates)

