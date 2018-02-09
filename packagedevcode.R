##############################################################################
# 
# package development code for LDATS
#
# Two Stage Analysis: Latent Dirichlet Allocation - Time Series 
#

  library(devtools)

  pkg_loc <- getwd()

  pkg_depends <- c("topicmodels", "RCurl", "multipanelfigure", "reshape2",  
                   "dplyr", "memoise", "lubridate", "progress", "ggplot2",  
                   "viridis", "nnet", "RColorBrewer", "Rcpp", "bindrcpp", 
                   "tidyverse", "gridExtra", "topicmodels", "doParallel")
  n_pkgs <- length(pkg_depends)

  for(i in 1:n_pkgs){
    devtools::use_package(pkg_depends[i], "Imports", pkg_loc)
  }

  devtools::load_all(devtools::as.package(pkg_loc))

  devtools::document(devtools::as.package(pkg_loc))

  # special definitions

    `%>%` <- dplyr::`%>%`
    `period` <- lubridate::`period`
    `%dopar%` <- foreach::`%dopar%`
