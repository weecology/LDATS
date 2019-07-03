library(testthat)
library(LDATS)
test_check("LDATS")

# note: if the figure images need to be updated, delete the files, set tenv to
#  thing other than "cran" in the relevant plot testing scripts, then run
#  vdiffr::manage_cases() and validate everything 
