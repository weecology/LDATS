# Vignettes that depend on internet access have been precompiled:

knitr::knit("vignettes/paper-comparison.Rmd.remote_needs", output = "vignettes/paper-comparison.Rmd")

# Manually move figure files from /figure to /vignettes/figure

dir.create(file.path("vignettes", "figure"), showWarnings = FALSE)
file.copy(file.path("figure"), file.path("vignettes"), recursive = TRUE)
unlink(file.path("figure"), force = TRUE, recursive = TRUE)

# Build checking

devtools::build_vignettes( )

