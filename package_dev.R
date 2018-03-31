pkg_depends <- c("topicmodels", "RCurl", "multipanelfigure", "reshape2",  
                 "dplyr", "memoise", "lubridate", "progress", "ggplot2",  
                 "viridis", "nnet", "RColorBrewer", "Rcpp", "bindrcpp", 
                 "tidyverse", "gridExtra", "topicmodels", "doParallel",
                 "doRNG", "rlang")
n_pkgs <- length(pkg_depends)
for(i in 1:n_pkgs){
  devtools::use_package(pkg_depends[i], "Imports")
}

devtools::load_all()
devtools::document()
#devtools::build()




