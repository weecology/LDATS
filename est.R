
package_TSx <- function(data, formula, timename, weights, control, rho_dist, 
                         eta_dist){

  if (is.null(nchangepoints)){
    nchangepoints <- 0
    mod <- multinom_TS(data, formula, changepoints = NULL, timename, weights,
                       control)
    mod <- mod[[1]][[1]]
    lls <- as.numeric(logLik(mod))
    rhos <- NULL
  } else{
    lls <- rho_dist$lls[1, ]
    rhos <- t(array(rho_dist$cpts[ , 1, ], dim = dim(rho_dist$cpts)[c(1, 3)]))
  }



  logLik <- mean(lls)
  ncoefs <- ncol(eta_dist)
  nparams <- nchangepoints + ncoefs 
  AIC <- -2 * logLik + 2 * nparams

  out <- list(lls = lls, rhos = rhos,
              etas = eta_dist, ptMCMC_diagnostics = ptMCMC_diagnostics,
              rho_summary = rho_summary, rho_vcov = rho_vcov,
              eta_summary = eta_summary, eta_vcov = eta_vcov,
              logLik = logLik, nparams = nparams, AIC = AIC)





}
