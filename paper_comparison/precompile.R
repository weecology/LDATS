# Vignettes that depend on internet access have been precompiled:

wd <- getwd()
setwd("vignettes/")

knitr::knit("paper-comparison.Rmd.remote_needs", output = "paper-comparison.Rmd")
setwd(wd)

# Build checking

devtools::build_vignettes( )

