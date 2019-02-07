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

## ----lda_set, eval =F----------------------------------------------------
#  lda_model_set <- LDA_set(document_term_table = rodents$document_term_table,
#                           topics = c(2:6),
#                           nseeds = 100,
#                          control = LDA_controls_list(quiet = TRUE))
#  

## ----lda set not quiet---------------------------------------------------
lda_model_set2 <- LDA_set(document_term_table = rodents$document_term_table,
                         topics = c(2:3),
                         nseeds = 10)

