devtools::load_all()
data(rodents)


x<-(c(rnorm(1e4, 0,1), rnorm(1e4, 10, 2)))
x
density(x)




document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
LDA_models <- LDA_set(document_term_table, topics = 2)[[1]]
data <- document_covariate_table
data$gamma <- LDA_models@gamma
weights <- document_weights(document_term_table)


TSmod <- TS(data, formula = gamma ~ 1, nchangepoints = 1, 
            timename = "newmoon", weights, control = list())



oh yeah and "multinom" now needs to be "compositional" :(
or how exactly?
categorical?


prob(x)^(1/T)
(1/T)*log(prob(x))



plot(exp((1*dnorm(seq(-1,1,0.01), 0, 1, log=TRUE))), ylim=c(0,0.4))
points(exp((20*dnorm(seq(-1,1,0.01), 0, 1, log=TRUE))))


plot(((1*dnorm(seq(-1,1,0.01), 0, 1, log=TRUE))), ylim=c(-4,0.4))
points(((4*dnorm(seq(-1,1,0.01), 0, 1, log=TRUE))))



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