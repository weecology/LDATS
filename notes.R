VCV
 what about fitting the lm models using nnet.default but with softmax = F
	oh whoa you get a block diag hessian...because uncor responses??
		yeah, looks like this is legit and ok
		when you run stuff without weights, you get nearly identical vcovs
		but the off diag 0s aren't truly 0
      yeah this is all v interesting, even in lm things align but with more
        actual noise among the response variables (still block diag tho)
	think that's likely the place to go for now
 will using the Hessian in the regular situation avoid the mirror issue?
predict
tidying check functions and general usage
examples
tests
tempering
vignettes
  

devtools::load_all()
   data(rodents)
LDAs <- LDA(data = rodents, topics = 3, replicates = c(1))

Y <- (LDAs[[1]][[1]]$document_topic_table)
Y2 <- ilr(Y)
X <- rodents$document_covariate_table[,"newmoon",drop=F]
weights <- document_weights(rodents$document_term_table)

mod_wlm_ilr <- lm(Y2~X$newmoon, weights = weights)

round(vcov(mod_lm_ilr), 6)
round(vcov(summary(mod_wlm_ilr)[[1]]), 6)
round(vcov(summary(mod_wlm_ilr)[[2]]), 6)



TSs <- TS(LDAs = LDAs, formulas = ~ newmoon, nchangepoints = 0, 
          timename = "newmoon", weights = TRUE,
          control = list(response = simplex_TS, 
                         method_args = 
                         list(control = ldats_classic_control(nit = 100)),
                                response_args = 
                           list(control = 
                            simplex_TS_control(transformation = "ilr"))))






mod_multinom <- multinom(Y~X$newmoon, Hess = TRUE, weights = weights)
round(vcov(mod_wmultinom), 6)


X2 <- model.matrix(mod_multinom)

mask <- c(F, F, F, F, T, T, F, T, T)

mod_wnnet <- nnet::nnet.default(X2, Y, mask = mask, weights = weights,
                           size = 0, rang = 0, skip = TRUE, softmax = TRUE,
                           Hess = TRUE, censored = FALSE)
round(solve(mod_wnnet$Hessian[mask, mask]), 6)



mod_lm_ilr <- lm(Y2~X$newmoon)
round(vcov(mod_lm_ilr), 6)

mod_wlm_ilr <- lm(Y2~X$newmoon, weights = weights)
round(vcov(summary(mod_wlm_ilr)[[1]]), 6)
round(vcov(summary(mod_wlm_ilr)[[2]]), 6)


mod_wnnet_ilr <- nnet::nnet(Y2~newmoon, data=X, weights = weights,
                           size = 0, skip = TRUE, softmax = FALSE,
                           Hess = TRUE, linout = TRUE)
round(solve(mod_wnnet_ilr$Hessian), 7)






WW <- rep(0, 6)
mask2 <- c(F, T, T, F, T, T)
mod_nnetd_ilr <- nnet::nnet.default(X2, Y2, weights = weights, Wts = WW,
                           size = 0, skip = TRUE, softmax = FALSE,
                           Hess = TRUE, linout = TRUE)
round(solve(mod_nnetd_ilr$Hessian[mask2, mask2]), 7)






TSs <- TS(LDAs = LDAs, formulas = ~ newmoon, nchangepoints = 0:1, 
          timename = "newmoon", weights = TRUE,
          control = list(response = simplex_TS, 
                         method_args = 
                         list(control = ldats_classic_control(nit = 100)),
                                response_args = 
                           list(control = 
                            simplex_TS_control(transformation = "alr"))))


LDATSs <- LDA_TS(data = rodents, topics = 2:3, replicates = c(1, 4),
                 formulas = ~ 1, nchangepoints = 1, timename = "newmoon",
                 weights = TRUE,
                 control = list(response = simplex_TS,
                                TS_method_args = 
                           list(control = ldats_classic_control(nit = 100)),
                                response_args = 
                           list(control = 
                            simplex_TS_control(transformation = "alr"))))


plot(LDAs)
plot(TSs)
plot(TSs, plot_type="diagnostic")
plot(LDATSs)


not sure if needed:

time_order_data <- function(x, timename = "time"){
  time_order <- order(x[ , timename])
  x[time_order , ]
}