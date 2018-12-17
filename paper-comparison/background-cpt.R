library(LDATS)

load('/Users/renatadiaz/Desktop/paper_cpt_ldats_lda.Rdata')
source('/Users/renatadiaz/Documents/GitHub/weecology/LDATS/paper-comparison/christensen-ecology-files/changepointmodel.r')
# set up parameters for model
year_continuous = 1970 + as.integer(julian(ch_dates)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)

# run models with 1, 2, 3, 4, 5 changepoints
cp_results_rodent = changepoint_model(ldats_lda_selected[[1]], x, 1, weights = rep(1,length(year_continuous)))
cp_results_rodent2 = changepoint_model(ldats_lda_selected[[1]], x, 2, weights = rep(1,length(year_continuous)))
cp_results_rodent3 = changepoint_model(ldats_lda_selected[[1]], x, 3, weights = rep(1,length(year_continuous)))
cp_results_rodent4 = changepoint_model(ldats_lda_selected[[1]], x, 4, weights = rep(1,length(year_continuous)))
cp_results_rodent5 = changepoint_model(ldats_lda_selected[[1]], x, 5, weights = rep(1,length(year_continuous)))
cp_results_rodent6 = changepoint_model(ldats_lda_selected[[1]], x, 6, weights = rep(1,length(year_continuous)))

save.image(file = '/Users/renatadiaz/Desktop/paper_cpt_ldats_lda.Rdata')