# working through an example based on the Portal rodent data

"%>%" <- dplyr::"%>%"
"period" <- lubridate::"period"
"%dopar%" <- foreach::"%dopar%"
"%dorng%" <- doRNG::"%dorng%"

r_data <- LDATS::create_rodent_table() 
lda_data <- r_data %>%
            dplyr::select(-c(period, censusdate, nplots, ntraps))
r_LDA <- LDATS::LDA(data = lda_data, ntopics = 2:5, nseeds = 2, ncores = 8)

ldamodel <- r_LDA[[1]]
ts_data <- data.frame(r_data[ , "censusdate"])
colnames(ts_data) <- "date"
ts_data[, 1] <- as.Date(ts_data[, 1])
diffyrs <- difftime(ts_data[ , 1], ts_data[1, 1], unit = "days") / 365.25
ts_data$diffyrs <- round(as.numeric(diffyrs), 2)
ts_data$gamma <- ldamodel@gamma

preds <- "1"
preds2 <- "diffyrs"
sample_sizes <- apply(lda_data, 1, sum)
wts <- sample_sizes/max(sample_sizes)
wts2 <- rep(1, length(wts))
cp <- NULL
cp2 <- as.Date("2000-01-01")
cp3 <- as.Date(c("1980-01-01", "2000-01-01"))

LDATS::mutilnom_ts(preds2, ts_data, changepoints = cp2, wts)




