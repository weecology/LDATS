---
title: "Comparison to Christensen et al. 2018"
author: "Renata Diaz"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{paper-comparison}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Side-by-side comparison of LDATS results with analysis from Christensen et al 2018. 

## LDATS Installation

To obtain the most recent version of **LDATS**, install the most recent 
version from GitHub:

```{r, eval=FALSE}
install.packages("devtools")
devtools::install_github("weecology/LDATS")
```

```{r, eval = T}
library(LDATS)

```

## Christensen 2018 analysis files

Download Christensen 2018 analysis scripts & data files from [Extreme-events-LDA repo:](https://github.com/emchristensen/Extreme-events-LDA)

Main Analysis Scripts:

  * rodent_LDA_analysis.R main script for analyzing rodent community change using LDA
  
  * rodent_data_for_LDA.R contains a function that creates the rodent data table used in analyses
  
  * AIC_model_selection.R contains functions for calculating AIC for different candidate LDA models
  
  * changepointmodel.r contains change-point model code
  
  * LDA-distance.R function for computing Hellinger distance analyses
    
Data:

  * Rodent_table_dat.csv table of rodent data, created by rodent_data_for_LDA.R

Figure scripts:

  * LDA_figure_scripts.R contains functions for making main plots in manuscript (Fig 1). Called from rodent_LDA_analysis.R


```{r download Christensen}

paper_filepath <- here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files')
test_filepath <- here::here('vignettes', 'christensenetal-ristensenetal-comparison-files', 'christensenetal-2018-files', 'rodent_LDA_analysis.r')

if(!file.exists(test_filepath)) {

files_to_download <- c('rodent_LDA_analysis.r', 'rodent_data_for_LDA.r', 'AIC_model_selection.R', 'changepointmodel.r', 'LDA-distance.R', 'Rodent_table_dat.csv', 'LDA_figure_scripts.R')

for(i in 1:length(files_to_download)) {
  download.file(url = paste0("https://raw.githubusercontent.com/emchristensen/Extreme-events-LDA/master/", files_to_download[i]),
                destfile = paste0(paper_filepath, '/', files_to_download[i]))
}

rm(files_to_download)
rm(i)

}

rm(test_filepath)
rm(paper_filepath)
```



## Data

The Portal rodents control data is included in the LDATS package:

```{r LDATS data}

data(rodents)

head(rodents[[1]])

```

Load the data used in Christensen et al: 

```{r Paper data}
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'rodent_data_for_LDA.r'))

dat = create_rodent_table(period_first = 1,
                          period_last = 436,
                          selected_plots = c(2,4,8,11,12,14,17,22),
                          selected_species = c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'))

# dates to go with count data
moondat = read.csv(text=getURL("https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/moon_dates.csv"),stringsAsFactors = F)
moondat$date = as.Date(moondat$censusdate)

period_dates = filter(moondat,period %in% rownames(dat)) %>% select(period,date)
dates = period_dates$date

paper_dat <- dat

paper_dates <- dates

rm(list = c('dat', 'dates', 'period_dates',
            'create_rodent_table'))

```


Compare paper data to LDATS data:

```{r rodent data comparison} 

compare <- rodents[[1]] == paper_dat

length(which(rowSums(compare) < ncol(compare)))
```


There are 16 rows where the data included in LDATS differs from the paper data. This is because the LDATS data is not adjusted to account for trapping effort, but the paper data has divided all census counts by the actual number of plots trapped and multiplied by 8 to account for incompletely-trapped censuses. 

To confirm this, refer to lines 36-46 in `rodent_data_for_LDA.r`:

```
# retrieve data on number of plots trapped per month
  trap_table = read.csv('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv')
  trap_table_controls = filter(trap_table, plot %in% selected_plots)
  nplots_controls = aggregate(trap_table_controls$sampled,by=list(period = trap_table_controls$period),FUN=sum)
  
  # adjust species counts by number of plots trapped that month
  r_table_adjusted = as.data.frame.matrix(r_table)
  for (n in 1:436) {
    #divide by number of control plots actually trapped (should be 8) and multiply by 8 to estimate captures as if all plots were trapped
    r_table_adjusted[n,] = round(r_table_adjusted[n,]/nplots_controls$x[n]*8)
  }
  
    return(r_table_adjusted)

```  
 
Running the same procedure on the LDATS data, we modify 16 rows, and finish with a data frame that matches the one from the paper.

```{r adjust LDATS data after Christensen et al}

trap_table = read.csv('https://raw.githubusercontent.com/weecology/PortalData/master/Rodents/Portal_rodent_trapping.csv')
  trap_table_controls =dplyr::filter(trap_table, plot %in% c(2,4,8,11,12,14,17,22))
  nplots_controls = aggregate(trap_table_controls$sampled,by=list(period = trap_table_controls$period),FUN=sum)
  
   # adjust species counts by number of plots trapped that month
  ldats_rodents_adjusted = as.data.frame.matrix(rodents[[1]])
  for (n in 1:436) {
    #divide by number of control plots actually trapped (should be 8) and multiply by 8 to estimate captures as if all plots were trapped
    ldats_rodents_adjusted[n,] = round(ldats_rodents_adjusted[n,]/nplots_controls$x[n]*8)
  }

  
compare_raw <- rodents[[1]] == ldats_rodents_adjusted
length(which(rowSums(compare_raw) < ncol(compare_raw)))

compare_adjusted <- ldats_rodents_adjusted == paper_dat
length(which(rowSums(compare_adjusted) < ncol(compare_adjusted)))

```

```{r data adjustment cleanup, include = FALSE}
rm(list = c('compare', 'compare_adjusted', 'compare_raw', 'n', 'ldats_rodents_adjusted', 'nplots_controls', 'trap_table', 'trap_table_controls'))
```


Because the LDA procedure weights the information from documents (census periods) according to the number of words (rodents captured), we now believe it is most appropriate to run the LDA on _unadjusted_ trapping data, and we recommend that users of LDATS do so. However, to maintain consistency with Christensen et al 2018, we will proceed using the _adjusted_ rodent table in this vignette. 

```{r switch to adjusted rodents}
rodents[[1]] <- paper_dat
```


The LDATS rodent data comes with a `document_covariate_table`, which we will use later as the predictor variables for the changepoint models. In this table, time is expressed as new moon numbers. Later we will want to be able to interpret the results in terms of census dates. We will add a column to the `document_covariate_table` to convert new moon numbers to census dates. We will not reference this column in any of the formulas we pass to the changepoint models, so it will be ignored until we need it.

```{r add dates to covariate table}

head(rodents$document_covariate_table)

moondat <- dplyr::select(moondat, newmoonnumber, censusdate)
colnames(moondat) <- c('newmoon', 'censusdate')

new_cov_table <- dplyr::left_join(rodents$document_covariate_table, moondat, by = 'newmoon')

rodents$document_covariate_table <- new_cov_table

```

```{r cleanup dates, include = F}
rm(moondat)
rm(new_cov_table)
```

## Stage 1: LDA

While LDATS can run start-to-finish with `LDATS::LDA_TS`, here we will work through the process function-by-function to isolate differences with the paper. For an illustration of `LDA_TS`, see the `codebase` vignette.

```{r LDATS LDAs, eval = F}

ldats_ldas <- LDATS::LDA_set(document_term_table = rodents$document_term_table, topics = c(2:6), nseeds = 100)
ldats_lda_selected <- LDATS::select_LDA(LDA_models = ldats_ldas)

save(ldats_ldas, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldas.Rds'))
save(ldats_lda_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))

```

```{r rm LDATS LDAS, include = F, eval = F}
rm(list = c('ldats_ldas', 'ldats_lda_selected'))
```

Paper LDAS:

```{r create dat for paper lda, include = F}
dat = paper_dat
```

```{r paper LDAs, eval = F}
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'AIC_model_selection.R'))
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'LDA-distance.R'))

# Fit a bunch of LDA models with different seeds
# Only use even numbers for seeds because consecutive seeds give identical results
seeds = 2*seq(200)

# repeat LDA model fit and AIC calculation with a bunch of different seeds to test robustness of the analysis
best_ntopic = repeat_VEM(paper_dat,
                         seeds,
                         topic_min=2,
                         topic_max=6)
hist(best_ntopic$k,breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5,9.5),xlab='best # of topics', main='')

# 2b. how different is species composition of 4 community-types when LDA is run with different seeds?
# ==================================================================
# get the best 100 seeds where 4 topics was the best LDA model
seeds_4topics = best_ntopic %>% 
  filter(k == 4) %>% 
  arrange(aic) %>% 
  head(100) %>% 
  pull(SEED)

# choose seed with highest log likelihood for all following analyses
#    (also produces plot of community composition for 'best' run compared to 'worst')
best_seed = calculate_LDA_distance(paper_dat,seeds_4topics)
mean_dist = unlist(best_seed)[2]
max_dist = unlist(best_seed)[3]

# ==================================================================
# 3. run LDA model
# ==================================================================
ntopics = 4
SEED = unlist(best_seed)[1]  # For the paper, use seed 206
ldamodel = LDA(paper_dat,ntopics, control = list(seed = SEED),method='VEM')

save(ldamodel,file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))
```

```{r cleanup paper LDAS, include = F, eval = F}
rm(list = c('dat', 'ldamodel', 'max_dist', 'mean_dist',
            'ntopics', 'SEED', 'seeds', 'seeds_4topics', 'aic_model', 'aic_model_gibbs', 'calculate_LDA_distance', 'Hellinger', 'min_H', 'repeat_VEM', 'best_ntopic', 'best_seed'))
```

```{r reload LDAS, include = FALSE}

load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))

```

### Plot LDAS

#### Paper LDA
```{r plot paper LDA, fig.width=6, fig.height=6}
# Paper
plot(ldamodel, cols = NULL, option = "D")

```

#### LDATS LDA
```{r plot LDATS LDA, fig.width=6, fig.height=6}
# LDATS
plot(ldats_lda_selected[[1]], cols = NULL, option = "D")
```

```{r remove LDAS after plotting, include = F}
rm(ldamodel)
rm(ldats_lda_selected)
```

The paper method finds 4 topics and LDATS finds 6. This is because of an update to the model selection procedure. The paper conservatively overestimates the number of parameters and therefore overpenalizes AIC for models with more topics. For this vignette, we will compare the results from using both LDA models. 


## Changepoint models

We will compare four combinations of LDA + changepoint models:

* LDATS LDA + LDATS changepoint
* LDATS LDA + paper changepoint
* Paper LDA + LDATS changepoint
* Paper LDA + paper changepoint

The paper changepoint model weighted all sample periods equally, wheras LDATS can weight sample periods according to how many individuals were captured (controlled by the `document_term_weights`). We now believe it is more appropriate to weight periods proportional to captures, and this is what we recommend for LDATS users. For the purposes of comparison, we will continue set all weights = 1 for both changepoint models. For an example of LDATS run with proportional weights, see the rodents vignette. [?]

### Running paper changepoint models

#### LDATS LDA and paper changepoint

```{r load LDATS LDA for paper cpt}

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))

```

```{r run LDATS LDA and paper cpt, eval = F}
#### Run changepoint ####
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))

# set up parameters for model
year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# run models with 1, 2, 3, 4, 5, 6 changepoints
cp_results_rodent = changepoint_model(ldats_lda_selected[[1]], x, 1, weights = rep(1,length(year_continuous)))
save(cp_results_rodent, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_1.Rds'))
rm(cp_results_rodent)

cp_results_rodent2 = changepoint_model(ldats_lda_selected[[1]], x, 2, weights = rep(1,length(year_continuous)))
save(cp_results_rodent2, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_2.Rds'))
rm(cp_results_rodent2)

cp_results_rodent3 = changepoint_model(ldats_lda_selected[[1]], x, 3, weights = rep(1,length(year_continuous)))
save(cp_results_rodent3, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_3.Rds'))
rm(cp_results_rodent3)

cp_results_rodent4 = changepoint_model(ldats_lda_selected[[1]], x, 4, weights = rep(1,length(year_continuous)))
save(cp_results_rodent4, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
rm(cp_results_rodent4)

cp_results_rodent5 = changepoint_model(ldats_lda_selected[[1]], x, 5, weights = rep(1,length(year_continuous)))
save(cp_results_rodent5, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_5.Rds'))
rm(cp_results_rodent5)

cp_results_rodent6 = changepoint_model(ldats_lda_selected[[1]], x, 6, weights = rep(1,length(year_continuous)))
save(cp_results_rodent6, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_6.Rds'))
rm(cp_results_rodent6)
```

```{r select LDATS LDA and paper cpt}
#### Changepoint model selection ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_1.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_2.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_3.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_5.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_6.Rds'))

ntopics = ldats_lda_selected[[1]]@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(6+1)+(6))

# The lowest deviance is for 4 changepoints. 


df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>% reshape::melt()
year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
df_4$value = year_continuous[df_4$value]


ldats_paper_cpt_dates <- vector(length = 4)
ldats_paper_cpt_dates[1] <- mean(df_4[df_4$variable=='V1','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[2] <-mean(df_4[df_4$variable=='V2','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[3] <-mean(df_4[df_4$variable=='V3','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
ldats_paper_cpt_dates[4] <-mean(df_4[df_4$variable=='V4','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')

ldats_paper_cpt_dates

```

```{r clean up LDATS LDA and paper changepoint, include = F}
rm(list = c('cp_results_rodent', 'cp_results_rodent2', 'cp_results_rodent3', 'cp_results_rodent4', 'cp_results_rodent5', 'cp_results_rodent6', 'ntopics', 'ldats_lda_selected', 'df_4'))
```

#### Paper LDA and paper changepoint

```{r load paper LDA for paper cpt}

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))

```

```{r run paper LDA and paper cpt, eval = F}
#### Run changepoint ####
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))

# set up parameters for model
year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# run models with 1, 2, 3, 4, 5, 6 changepoints
cp_results_rodent = changepoint_model(ldamodel, x, 1, weights = rep(1,length(year_continuous)))
save(cp_results_rodent, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper1.Rds'))
rm(cp_results_rodent)

cp_results_rodent2 = changepoint_model(ldamodel, x, 2, weights = rep(1,length(year_continuous)))
save(cp_results_rodent2, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper2.Rds'))
rm(cp_results_rodent2)

cp_results_rodent3 = changepoint_model(ldamodel, x, 3, weights = rep(1,length(year_continuous)))
save(cp_results_rodent3, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper3.Rds'))
rm(cp_results_rodent3)

cp_results_rodent4 = changepoint_model(ldamodel, x, 4, weights = rep(1,length(year_continuous)))
save(cp_results_rodent4, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper4.Rds'))
rm(cp_results_rodent4)

cp_results_rodent5 = changepoint_model(ldamodel, x, 5, weights = rep(1,length(year_continuous)))
save(cp_results_rodent5, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper5.Rds'))
rm(cp_results_rodent5)

cp_results_rodent6 = changepoint_model(ldamodel, x, 6, weights = rep(1,length(year_continuous)))
save(cp_results_rodent6, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper6.Rds'))
rm(cp_results_rodent6)
```

```{r select paper LDA and paper cpt}
#### Changepoint model selection ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_1.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_2.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_3.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_4.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_5.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_6.Rds'))

ntopics = ldamodel@k

# change point model selection
# mean deviance ( -2 * log likelihood) + 2*(#parameters)
mean(cp_results_rodent$saved_lls * -2) + 2*(3*(ntopics-1)*(1+1)+(1))
mean(cp_results_rodent2$saved_lls * -2)+ 2*(3*(ntopics-1)*(2+1)+(2))
mean(cp_results_rodent3$saved_lls * -2)+ 2*(3*(ntopics-1)*(3+1)+(3))
mean(cp_results_rodent4$saved_lls * -2)+ 2*(3*(ntopics-1)*(4+1)+(4))
mean(cp_results_rodent5$saved_lls * -2)+ 2*(3*(ntopics-1)*(5+1)+(5))
mean(cp_results_rodent6$saved_lls * -2)+ 2*(3*(ntopics-1)*(6+1)+(6))

# The lowest deviance is for 4 changepoints. 

df_4 = as.data.frame(t(cp_results_rodent4$saved[,1,])) %>%  reshape::melt()
df_4$value = year_continuous[df_4$value]


paper_paper_cpt_dates <- vector(length = 4)
paper_paper_cpt_dates[1] <- mean(df_4[df_4$variable=='V1','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[2] <-mean(df_4[df_4$variable=='V2','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[3] <-mean(df_4[df_4$variable=='V3','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')
paper_paper_cpt_dates[4] <-mean(df_4[df_4$variable=='V4','value']) %>% lubridate::date_decimal() %>% format('%Y-%m-%d')

paper_paper_cpt_dates

```

```{r clean up paper LDA and paper changepoint, include = F}
rm(list = c('cp_results_rodent', 'cp_results_rodent2', 'cp_results_rodent3', 'cp_results_rodent4', 'cp_results_rodent5', 'cp_results_rodent6', 'ntopics', 'ldamodel', 'df_4'))
```

### Running LDATS changepoint models

#### LDATS LDA and LDATS changepoint

```{r load LDATS LDA for LDATS cpt}

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_lda_selected.Rds'))

```
#### Run the TS models

```{r run LDATS LDA and LDATS changepoint, eval = F}
#### Run LDATS changepoint ####

ldats_ldats_cpt <- TS_on_LDA(LDA_models = ldats_lda_selected, 
                             document_covariate_table = rodents$document_covariate_table,
                             formulas = ~ sin_year + cos_year,
                             nchangepoints = 1:6,
                             weights = NULL)


save(ldats_ldats_cpt, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats.Rds'))
rm(ldats_ldats_cpt)
rm(ldats_lda_selected)
```

```{r load LDATS LDA and LDATS changepoint, include = F}
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats.Rds'))
```

### Select model

```{r select ldats lda + ldats cpt}

ldats_ldats_cpt_selected <- select_TS(ldats_ldats_cpt)

ldats_ldats_cpt_selected$nchangepoints

```

### Get dates associated with changepoints

```{r get ldats lda + ldats cpt dates}
ldats_ldats_cpt_dates <- ldats_ldats_cpt_selected$rho_summary$Mean
ldats_ldats_cpt_dates <- round(ldats_ldats_cpt_dates)

ldats_ldats_cpt_dates <- as.data.frame(ldats_ldats_cpt_dates)
colnames(ldats_ldats_cpt_dates) <- 'newmoon'
ldats_ldats_cpt_dates <- dplyr::left_join(ldats_ldats_cpt_dates, rodents$document_covariate_table, by = 'newmoon')
ldats_ldats_cpt_dates <- ldats_ldats_cpt_dates$censusdate

```

```{r cleanup LDATS LDA and LDATS changepoint, include = F}

save(ldats_ldats_cpt_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats_selected.Rds'))

rm(ldats_ldats_cpt_selected)
rm(ldats_ldats_cpt)
```

### Paper LDA and LDATS changepoint


```{r load paper LDA for LDATS cpt}

#### Load LDA ####
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldamodel.Rds'))

```
#### Run the TS models

```{r run paper LDA and LDATS changepoint, eval = F}
#### Run LDATS changepoint ####

paper_ldats_cpt <- TS_on_LDA(LDA_models = ldamodel, 
                             document_covariate_table = rodents$document_covariate_table,
                             formulas = ~ sin_year + cos_year,
                             nchangepoints = 1:6,
                             weights = NULL)


save(paper_ldats_cpt, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats.Rds'))
rm(paper_ldats_cpt)
rm(ldamodel)
```

```{r load paper LDA and LDATS changepoint, include = F}
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats.Rds'))
```

### Select model

```{r select paper lda + ldats cpt}

paper_ldats_cpt_selected <- select_TS(paper_ldats_cpt)

paper_ldats_cpt_selected$nchangepoints

```

#### Get dates

Unlike the paper changepoint model, LDATS can recognize that sampling periods may not be equidistant, and can place changepoint estimates at new moons if they fall between nonconsecutive sampling periods. We can estimate the dates corresponding to those new moons, extrapolating from the census dates for adjacent census periods. 

```{r get paper lda + ldats cpt dates}
paper_ldats_cpt_dates <- paper_ldats_cpt_selected$rho_summary$Mean
paper_ldats_cpt_dates <- round(paper_ldats_cpt_dates)

paper_ldats_cpt_dates <- as.data.frame(paper_ldats_cpt_dates)
colnames(paper_ldats_cpt_dates) <- 'newmoon'
paper_ldats_cpt_dates <- dplyr::left_join(paper_ldats_cpt_dates, rodents$document_covariate_table, by = 'newmoon')

halfinterval = as.integer((as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[1] + 1)]) - as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[1] - 1)]))/2)

newdate <- lubridate::ymd(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[1] -1 )]) + lubridate::days(halfinterval)

paper_ldats_cpt_dates$censusdate[1] <- as.character(newdate)



halfinterval =as.integer((as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[3] + 2)]) - as.Date(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[3] - 1)]))/3)

newdate <- lubridate::ymd(rodents$document_covariate_table$censusdate[ which(rodents$document_covariate_table$newmoon == paper_ldats_cpt_dates$newmoon[3] -2 )]) + lubridate::days(halfinterval)

paper_ldats_cpt_dates$censusdate[3] <- as.character(newdate)

paper_ldats_cpt_dates <- paper_ldats_cpt_dates$censusdate
```


```{r cleanup paper LDA and LDATS changepoint, include = F}

save(paper_ldats_cpt_selected, file = here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats_selected.Rds'))

rm(paper_ldats_cpt_selected)
rm(paper_ldats_cpt)
rm(list = c('newdate', 'halfinterval'))
```

All of the models find four changepoints.


### Plot changepoint models

```{r load cpts to plot, include = F}
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_ldats_selected.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_ldats_selected.Rds'))
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'paper_paper_4.Rds'))
paper_paper_cpt_selected <- cp_results_rodent4
load(here::here('vignettes', 'christensenetal-comparison-files', 'model-cache', 'ldats_paper_4.Rds'))
ldats_paper_cpt_selected <- cp_results_rodent4

year_continuous = 1970 + as.integer(julian(paper_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)
ntopics = 4

```

#### Paper LDA and LDATS changepoint

```{r plot paper LDA and LDATS cpts, fig.width=6, fig.height=6}

plot(paper_ldats_cpt_selected)

```

#### LDATS LDA and LDATS changepoint
```{r plot ldats LDA and LDATS cpt, fig.width=6, fig.height=6}
plot(ldats_ldats_cpt_selected)
```

#### Paper LDA and paper changepoint
```{r plot paper LDA and paper cpt, fig.width=6, fig.height=6}

source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'LDA_figure_scripts.R'))
source(here::here('vignettes', 'christensenetal-comparison-files', 'christensenetal-2018-files', 'changepointmodel.r'))
paper_cpts = find_changepoint_location(paper_paper_cpt_selected)
paper_cpt_plot = get_ll_non_memoized_plot(ldamodel,x,paper_cpts,make_plot=T,weights=rep(1,length(year_continuous)))
annual_hist(paper_paper_cpt_selected, year_continuous)
paper_cpt_plot

```


#### LDATS LDA and paper changepoint
```{r plot LDATS lda and paper cpt, fig.width=6, fig.height=6}
ntopics = 6
ldats_cpts = find_changepoint_location(ldats_paper_cpt_selected)
ldats_cpt_plot = get_ll_non_memoized_plot(ldats_lda_selected[[1]],x,ldats_cpts,make_plot=T,weights=rep(1,length(year_continuous)))
annual_hist(ldats_paper_cpt_selected, year_continuous)
ldats_cpt_plot
```


The results of the changepoint model appear extremely robust to both choice of LDA model and choice of changepoint model. 

### Report changepoint dates
```{r report cpt dates, include = F}

cpt_dates <- as.data.frame(paper_paper_cpt_dates)
cpt_dates <- cbind(cpt_dates, ldats_paper_cpt_dates)
cpt_dates <- cbind(cpt_dates, ldats_ldats_cpt_dates)
cpt_dates <- cbind(cpt_dates, paper_ldats_cpt_dates)

colnames(cpt_dates) <- c('paperLDA_papercpt', 'ldatsLDA_papercpt', 'ldatsLDA_ldatscpt', 'paperLDA_ldatscpt')

```

```{r print cpt dates}
cpt_dates
```

The choice of LDA has more influence on the changepoint locations than the choice of changepoint model - probably because the LDATS LDA has 6 topics, and the paper LDA has 4. However, all of the models agree to within 6 months in most cases, and a year for the broader early 1990s changepoint. 

```{r cleanup cpt_dates, include = F}
rm(list= 'ldats_ldats_cpt_dates', 'ldats_paper_cpt_dates',  'paper_ldats_cpt_dates', 'paper_paper_cpt_dates', 'cpt_dates')
```
