##############################################################################
# 
# package development code for LDATS
#
# Two Stage Analysis: Latent Dirichlet Allocation - Time Series 
#
# held under MIT
#

  library(devtools)

  pkg_loc <- getwd()

  pkg_depends <- c("topicmodels", "RCurl", "multipanelfigure", "reshape2",  
                   "dplyr", "memoise", "lubridate", "progress", "ggplot2",  
                   "viridis", "nnet", "RColorBrewer", "Rcpp", "bindrcpp", 
                   "tidyverse", "gridExtra", "topicmodels", "doParallel")

  for(i in 1:length(pkgdpns)){
    devtools::use_package(pkg_depends[i], "Imports", pkgloc)
  }

  devtools::load_all(devtools::as.package(pkgloc))

  devtools::document(devtools::as.package(pkgloc))

  # define the pipe operator, period function, and dopar operator 

    `%>%` <- dplyr::`%>%`
    `period` <- lubridate::`period`
    `%dopar%` <- foreach::`%dopar%`
