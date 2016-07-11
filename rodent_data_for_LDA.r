# code used for getting rodent data into table form for LDA analysis



rod = read.csv('../PortalData/Rodents/Portal_rodent.csv',na.strings = '',as.is=T)

# remove negative period, non-census
rod = rod[rod$period>0,]

# only up to period 436, when plots were switched
rod = rod[rod$period<=436,]

# controls only
controls = c(2,4,8,11,12,14,17,22)
rod = rod[rod$plot %in% controls,]

# target species only
rod = rod[rod$species %in% c('BA','DM','DO','DS','NA','OL','OT','PB','PE','PF','PH','PI','PL','PM','PP','RF','RM','RO','SF','SH','SO'),]

# aggregate by period and species
r_table = table(rod$period,rod$species)

# adjust for number of plots trapped that month
trap_table = read.csv('../PortalData/Rodents/Portal_rodent_trapping.csv')
trap_table_controls = trap_table[trap_table$Plot %in% controls,]
nplots_controls = aggregate(trap_table_controls$Sampled,by=list(Period = trap_table_controls$Period),FUN=sum)

r_table_adjusted = as.data.frame.matrix(r_table)
for (n in 1:436) {
  #divide by number of control plots actually trapped (should be 8) and multiply by 8 to estimate captures as if all plots were trapped
  r_table_adjusted[n,] = round(r_table_adjusted[n,]/nplots_controls$x[n]*8)
}

write.table(r_table_adjusted,'Rodent_table_dat.csv',sep=',',row.names = F)
