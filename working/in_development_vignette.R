data(rodents)
?rodents
library(dplyr)
lda_data <- rodents %>%
            select(-c(newmoon, date, plots, traps))
ts_data <- rodents %>%
           select(c(newmoon)) %>% 
           rename(time = newmoon)
          

weights <- doc_weights(lda_data)

ldas <- parLDA(lda_data, ntopics = 2:3) 
selected <- LDA_select(ldas) 

prepped <- selected %>% 
           MTS_prep(ts_data)


mtss <- selected %>% 
        LDATS::MTS_prep(ts_data) %>%
        LDATS::MTS_set(formula = ~1, nchangepoints = 1, weights) 

xx<-LDATS::MTS(data[[1]], formula, nchangepoints, weights, nit = 10)