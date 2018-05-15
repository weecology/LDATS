data(rodents)
?rodents

lda_data <- rodents %>%
            dplyr::select(-c(newmoon, date, plots, traps))
ts_data <- rodents %>%
           dplyr::select(c(newmoon)) %>% 
           dplyr::rename(time = newmoon)
          

weights <- LDATS::doc_weights(lda_data)

ldas <- LDATS::LDA(lda_data, ntopics = 2:3) 
selected <- LDATS::LDA_select(ldas) 

prepped <- selected %>% 
        LDATS::MTS_prep(ts_data)


mtss <- selected %>% 
        LDATS::MTS_prep(ts_data) %>%
        LDATS::MTS_set(formula = ~1, nchangepoints = 1, weights) 

xx<-LDATS::MTS(data[[1]], formula, nchangepoints, weights, nit = 10)