# Do the data tables match?

ch_dat <- read.csv('rmd-ignore/previous-work/paper_dat.csv', 
                    stringsAsFactors = F)

library(LDATS)
data(rodents)
ldats_dat <- rodents$document_term_table

compare <- ldats_dat == ch_dat

which(!compare)
ldats_dat[351,]
ch_dat[351, ]
