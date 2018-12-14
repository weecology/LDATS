# to do:
#   tidy up documentation
#   full vignette
#
# working through documentation right now
#  scripts to do
#   ptMCMC
#   TS
#   TS_on_LDA
#   TS_plots
#   utilities

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


\ifelse{html}{\out{<b><i></i></b>}}{\eqn{\mathbf{}}}
\ifelse{html}{\out{}}{\eqn{\boldsymbol{}}}

#'   (\ifelse{html}{\out{<i>k<sub>m<sub>1</sub></sub></i>}}{\eqn{k_{m_1}}} 
#'   (\ifelse{html}{\out{<i>P<sub>m<sub>2</sub></sub></i>}}{\eqn{P_{m_2}}} 

#'   \ifelse{html}{
#'      \out{<i>&rho;<sub>m<sub>2</sub></sub></i>}}{
#'      \eqn{\boldsymbol{\rho}_{m_2}}} 
#'


#'   \ifelse{html}{\out{<b><i>t</i></b>}}{\eqn{\mathbf{t}}} 

#'   \ifelse{html}{\out{<b><i>v</i></b>}}{\eqn{\mathbf{v}}} 

#'   \ifelse{html}{\out{<b><i>t</i></b>}}{\eqn{\mathbf{t}}}
#'   \ifelse{html}{\out{<b><i>X</i></b>}}{\eqn{\mathbf{X}}}

#'   \ifelse{html}{
#'      \out{<span style="text-decoration: overline"><b><i>&Gamma;
#'           </i></b></span>}}{\eqn{\overline{\boldsymbol{\Gamma}}}}

#'   \ifelse{html}{
#'      \out{<span style="text-decoration: overbrace"><b><i>&rho;
#'           </i></b></span>}}{\eqn{\overbrace{\boldsymbol{\rho}}}}