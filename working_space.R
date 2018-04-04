# working through an example based on the Portal rodent data

"%>%" <- dplyr::"%>%"
"period" <- lubridate::"period"
"%dopar%" <- foreach::"%dopar%"
"%dorng%" <- doRNG::"%dorng%"

r_data <- LDATS::create_rodent_table() 
lda_data <- r_data %>%
            dplyr::select(-c(newmoonnumber, newmoondate, nplots, ntraps))
r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 8)

ldamodel <- r_LDA[[5]]

newmoondate <- as.Date(r_data$newmoondate)
nye <- as.Date(paste(format(newmoondate, "%Y"), "-12-31", sep = ""))
jday <- as.numeric(format(newmoondate, "%j"))
nye_jday <- as.numeric(format(nye, "%j"))
fr_of_yr <- round(jday / nye_jday, 3)

ts_data <- data.frame(time = r_data$newmoonnumber, foy = fr_of_yr)
ts_data$gamma <- round(ldamodel@gamma, 3)

preds <- "1"
preds2 <- "time + I(time^2)"
sample_sizes <- apply(lda_data, 1, sum)
wts <- round(sample_sizes/max(sample_sizes), 3)

LDATS::MTS(formula = preds, data = ts_data, nchangepoints = 1, weights = wts, 
  nit = 100)
