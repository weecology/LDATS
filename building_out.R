
system.time(for(i in 1:1e2){

    selection <- cbind(pdist$which_kicked[i, ], 1:ntemps)
    prop_changepts <- changepts
    curr_changepts_s <- changepts[selection]
    prop_changepts_s <- curr_changepts_s + pdist$kicks[i, ]
    #prop_changepts[selection] <- prop_changepts_s #HERE
    #prop_lls <- lls

    for (j in 1:ntemps){
      mod <- multinom_ts(formula, data, prop_changepts[ , j], weights)
      prop_lls[j] <- tryCatch(mod$logLik, error = function(x) {-Inf})
    }
    accepts <- runif(ntemps) < exp((prop_lls - lls) * betas)
    accept_rate <- accept_rate + accepts / nit
    changepts[ , accepts] <- prop_changepts[ , accepts]
    lls[accepts] <- prop_lls[accepts]
    
})




system.time(for(W in 1:1e3){



    for (j in 1:ntemps){
      mod <- ts_memo(formula, data, prop_changepts[ , j], weights)
      prop_lls[j] <- mod$logLik
    }
    accepts <- runif(ntemps) < exp((prop_lls - lls) * betas)
    accept_rate <- accept_rate + accepts / nit
    changepts[ , accepts] <- prop_changepts[ , accepts]
    lls[accepts] <- prop_lls[accepts]
    
})