

devtools::load_all()
data(rodents)
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table


xx <- LDA_TS(document_term_table, document_covariate_table,
                   topics = 2, nseeds = 1, formulas = ~ 1, nchangepoints = 1,
                   weights = document_weights(document_term_table), 
                   control = LDA_TS_controls_list())

xx <- LDA_TS(document_term_table, document_covariate_table,
                   topics = 4, nseeds = 10, formulas = ~ sin_year + cos_year, 
                   nchangepoints = 4,
                   weights = document_weights(document_term_table), 
                   control = LDA_TS_controls_list(
                             TS_control = TS_controls_list(nit = 1e3)))

#'   This function was designed to work within \code{\link{TS}} and process
#'   the output of \code{\link{est_changepts}}, but has been generalized
#'   and would work with any output from a ptMCMC as long as \code{ptMCMCout}
#'   is formatted properly.

#'   This function was designed to work within \code{\link{TS}} and 
#'   specifically \code{\link{est_changepts}}. It is still hardcoded to do
#'   so, but has the capacity to be generalized to work with any estimation
#'   via ptMCMC with additional coding work.
