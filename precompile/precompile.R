# Vignettes that depend on internet access are precompiled following
#    https://ropensci.org/blog/2019/12/08/precompute-vignettes/.

# temporarily change working directory 

wd <- getwd()
setwd("paper_comparison/")

knitr::knit(input  = "paper-comparison.Rmd.remote_needs", 
            output = "../vignettes/paper-comparison.Rmd")


# Move the output figure files to the vignettes folder

from_files <- list.files(path       = "output", 
                         full.names = TRUE)[grep(".png", list.files("output"))]

file.copy(from      = from_files,
          to        = file.path("..", "vignettes", "output"),
          recursive = TRUE,
          overwrite = TRUE)

file.copy(from      = file.path("figure"),
          to        = file.path("..", "vignettes"),
          recursive = TRUE,
          overwrite = TRUE)


# revert the working directory

setwd(wd)


# Build checking

# devtools::build_vignettes( )

