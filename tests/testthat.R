library(testthat)
library(LDATS)
tryCatch(dev.off(), error = function(x){})
test_check("LDATS")

# note: if the figure images need to be updated, delete the files, then run
#  vdiffr::manage_cases() and validate everything 