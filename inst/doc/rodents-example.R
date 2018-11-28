## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include=FALSE------------------------------------------------------
library(LDATS)
vers <- packageVersion("LDATS")
today <- Sys.Date()

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("weecology/LDATS")

## ------------------------------------------------------------------------
data(rodents)
head(rodents$document_term_table, 10)
head(rodents$document_covariate_table, 10)

