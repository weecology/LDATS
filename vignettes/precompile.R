# Vignettes that depend on internet access have been precompiled:

knitr::knit("vignettes/paper-comparison.Rmd.remote_needs", output = "vignettes/paper-comparison.Rmd")

# Manually move figure files from /figure to /vignettes


file.copy(from = list.files(file.path("figure"), full.names = TRUE),
          to   = file.path("vignettes"))

unlink(file.path("figure"), force = TRUE, recursive = TRUE)

# Build checking

devtools::build_vignettes( )

