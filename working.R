# have a rough set up for the data splitting, now need to actually 
# connect the pieces under the hood...like measurer? and selector?






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

# blech this destroys our ability to memoise tho
# i think we needed to cut that tie sooner or later tho tbh
# yeah we'll need to set up some of the interface at this level to be 
# much more flexed 
# create a function that acts like multinom for a basis-based model?
# or do we drop multinom all together and use alr for that?