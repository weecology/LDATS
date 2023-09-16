# Vignettes that depend on internet access have been precompiled:

wd <- getwd()
setwd(file.path(wd, "vignettes"))
knitr::knit("paper-comparison.Rmd.remote_needs", output = "paper-comparison.Rmd")


# Move figure files from /figure to /vignettes


file.copy(from = list.files(file.path("figure"), full.names = TRUE),
          to   = file.path("."))

unlink(file.path("figure"), force = TRUE, recursive = TRUE)

# Build checking

devtools::build_vignettes( )

