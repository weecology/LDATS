##############################################################################
# 
# package development code for LDA.pointbreak
#
# Latent Dirichlet Allocation with Break Point Analysis
#
# version 0.0.1 January 2018
#
# held under MIT
#
# #RequisiteJohnnyUtahReference
#
# Erica Christensen, David Harris, Hao Ye, and Juniper Simonis
#
##############################################################################

  # load devtools

    library(devtools)

  # set the package location

    pkgloc <- getwd()

  # add dependencies to the description and load them here 

    pkgdpns <- c("topicmodels", "RCurl", "multipanelfigure", "reshape2",  
                 "dplyr", "memoise", "lubridate", "progress", "ggplot2",  
                 "viridis", "nnet", "RColorBrewer", "Rcpp", "bindrcpp", 
                 "tidyverse", "gridExtra", "topicmodels", "doParallel")

    for(i in 1:length(pkgdpns)){
      devtools::use_package(pkgdpns[i], "Imports", pkgloc)
    }

  # load the package

    devtools::load_all(devtools::as.package(pkgloc))

  # populate the documentation

    devtools::document(devtools::as.package(pkgloc))

  # define the pipe operator, period function, and dopar operator 

    `%>%` <- dplyr::`%>%`
    `period` <- lubridate::`period`
    `%dopar%` <- foreach::`%dopar%`
