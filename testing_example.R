devtools::load_all()
data(rodents)
length(rodents)

rodents[["metadata"]] <- list(timename = "newmoon")
names(rodents) <- c("abundance", "covariates", "metadata")

rodents[["metadata"]] <- NULL
names(rodents)[1] <- "term"
mod <- LDA(rodents)
mod <- LDA(rodents, topics = 1)
mod <- LDA(rodents, control = LDA_control(LDA_args = list(method = "Gibbs")))

mod <- LDA(rodents[[1]])

