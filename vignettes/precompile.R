#Vignettes that depend on internet access have been precompiled:


knitr::knit("vignettes/paper-comparison.Rmd.remote_needs", output = "vignettes/paper-comparison.Rmd")

devtools::build_vignettes( )