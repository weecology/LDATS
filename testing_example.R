devtools::load_all()
data(rodents)
length(rodents)
rodents[["metadata"]] <- list(timename = "newmoon")

rodents[["metadata"]] <- NULL
mod <- LDA(rodents)
mod <- LDA(rodents, topics = 1)
mod <- LDA(rodents, control = LDA_control(LDA_args = list(method = "Gibbs")))