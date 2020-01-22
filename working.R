# need to remove restrictions on the selection functions:
#  if selector is NULL, return the set of models
#  if >1 model is selected, return all

devtools::load_all()
data(rodents)
names(rodents)



mod <- LDA_TS(data = rodents, topics = 2, nseeds = 1, formulas = ~1,
                 nchangepoints = 0, timename = "newmoon")

mmod <- meta_LDA_TS(data = rodents, topics = 2, nseeds = 1, formulas = ~1,
                 nchangepoints = 0, timename = "newmoon")

