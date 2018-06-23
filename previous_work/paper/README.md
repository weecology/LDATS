# Extreme-events-LDA
## Repo for holding code for working on LDA analysis of rodent data

This repo contains the analyses for my project on rodent community change, currently up as a preprint: https://doi.org/10.1101/163931

## Main Analysis Scripts: 
  * __rodent_LDA_analysis.R__ main script for analyzing rodent community change using LDA
  * __rodent_data_for_LDA.R__ contains a function that creates the rodent data table used in analyses
  * __AIC_model_selection.R__ contains functions for calculating AIC for different candidate LDA models
  * __changepointmodel.r__ contains change-point model code
  * __create_sim_data.R__ contains functions for creating simulated data, used in the manuscript supplement to demonstrate methods
  * __LDA-distance.R__ function for computing Hellinger distance analyses
  * __simulated_data_LDA.R__ LDA/changepoint analysis using simulated data, used in supplement to demonstrate methods
  
## Figure Scripts:
  * __NDVI_figure.R__ creates a timeseries of NDVI at the site from 1984-2015, to show periods of low resources (supplemental figure in manuscript)
  * __abundance_plots.R__ contains code for plotting timeseries of total rodent abundance (Fig 2 in manuscript)
  * __explanatory-figure.R__ code to create graphical representation of LDA model (supplemental fig in manuscript)
  * __LDA_figure_scripts.R__ contains functions for making main plots in manuscript (Fig 1). Called from rodent_LDA_analysis.R
  * __pie_plots_for_presentation.R__ creates pie plot of species composition at a given trapping time. Not used in manuscript, but I use it in various talks I've given
  
## Data files:
  * _Rodent_table_dat.csv_ table of rodent data, created by rodent_data_for_LDA.R
  * _Monthly_Landsat_NDVI.csv_ NDVI data for portal from Landsat satellites

## Folders:
  * __preliminaryscripts__ includes code from the initial stages of the project that didn't end up getting used in the final product
      * _abundance_plots_old.R_ is an old version of abundance_plots.R
      * _gibbs_functions.R_ contains functions for doing the gibbs sampler version of LDA (as opposed to VEM method)
      * _LDA_analysis.R_ early, incomplete version of rodent_LDA_analysis.R
      * _LDA_confidence_intervals.R_ early version of LDA analysis, playing with gibbs sampler and credible intervals
      * _rodent_abundance_histograms.R_ playing around with histograms of total rodent abundance, to show lowest abundance events
      * _rodent_euclid_dist_plot.R_ plots rodent community change over time, paired with a picture of shrub growth at Portal. This figure was in the first version of the manuscript but not later versions
      * _Thibault_2004_eucliddist.R_ playing around with the euclidean distance analysis from Thibault et al 2004
