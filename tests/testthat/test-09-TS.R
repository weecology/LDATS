context("Check TS functions")

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
topics <- 2
nseeds <- 1
formulas <- ~ 1
nchangepoints <- 1
weights <- document_weights(document_term_table)
LDAs <- LDA_set(document_term_table, topics, nseeds)
LDA_models <- select_LDA(LDAs)
control <- list(nit = 20, seed = 1)
timename <- "newmoon"
mods <- expand_TS(LDA_models, formulas, nchangepoints)
formula <- mods$formula[[1]]
nchangepoints <- mods$nchangepoints[1]
data <- prep_TS_data(document_covariate_table, LDA_models, mods, 1)
TSmod <- TS(data, formula, nchangepoints, timename, weights, control)

rho_dist <- est_changepoints(data, formula, nchangepoints, timename, weights,
                             control)
eta_dist <- est_regressors(rho_dist, data, formula, timename, weights,
                           control)
set.seed(1)
test_that("check measure_eta_vcov", {
  expect_is(measure_eta_vcov(eta_dist), "matrix")
  expect_equal(round(measure_eta_vcov(eta_dist)[1, 1], 2), 0.32)
  expect_error(measure_eta_vcov("ok"))
})

test_that("check measure_rho_vcov", {
  set.seed(1)
  nchangepoints <- dim(rho_dist$cpts)[1]
  if (is.null(nchangepoints)){
    nchangepoints <- 0
    mod <- multinom_TS(data, formula, changepoints = NULL, 
                       timename, weights, control)
    mod <- mod[[1]][[1]]
    lls <- as.numeric(logLik(mod))
    rhos <- NULL
  } else{
    lls <- rho_dist$lls[1, ]
    rhos <- t(array(rho_dist$cpts[ , 1, ], dim = dim(rho_dist$cpts)[c(1, 3)]))
  }
  expect_is(measure_rho_vcov(rhos), "matrix")
  expect_equal(round(measure_rho_vcov(rhos)[1, 1], 1), 572.4)
  expect_error(measure_rho_vcov("ok"))
})


test_that("check summarize_etas", {
  sum_e <- summarize_etas(eta_dist)
  expect_is(sum_e, "data.frame")
  expect_equal(round(sum_e[1, 1], 2), 2.52)
  expect_equal(summarize_etas(eta_dist[1:3, ])$AC10[1], as.factor("-"))
  expect_error(summarize_etas("ok"))
  expect_error(summarize_etas(eta_dist, LDA_controls_list()))
})

test_that("check summarize_rhos", {
  set.seed(1)
  nchangepoints <- dim(rho_dist$cpts)[1]
  if (is.null(nchangepoints)){
    nchangepoints <- 0
    mod <- multinom_TS(data, formula, changepoints = NULL, 
                       timename, weights, control)
    mod <- mod[[1]][[1]]
    lls <- as.numeric(logLik(mod))
    rhos <- NULL
  } else{
    lls <- rho_dist$lls[1, ]
    rhos <- t(array(rho_dist$cpts[ , 1, ], dim = dim(rho_dist$cpts)[c(1, 3)]))
  }

  sum_r <- summarize_rhos(rhos)
  expect_is(sum_r, "data.frame")
  expect_equal(round(sum_r[1, 1], 2), 253)
  expect_error(summarize_rhos("ok"))
  expect_error(summarize_rhos(rhos, LDA_controls_list()))
})



test_that("check est_changepoints", {
  set.seed(1)
  rhos <- est_changepoints(data, formula, nchangepoints, timename, weights,
                           control)
  expect_is(rhos, "list")
  expect_equal(length(rhos), 5)
  expect_equal(names(rhos), 
               c("cpts", "lls", "ids", "step_accepts", "swap_accepts"))
  expect_equal(dim(rhos[[1]]), c(1, 6, 20))
  expect_equal(round(sum(rhos$lls), 1), -29033.3)

  expect_error(est_changepoints(data, formula, nchangepoints,  
               timename, weights,"ok"))
  expect_error(est_changepoints(data, formula, nchangepoints, "ok", weights,
               control))
  expect_error(est_changepoints(data, formula, nchangepoints, "timename", 
               "ok", control))
  expect_error(est_changepoints(data, formula, "ok",   
               timename, weights, control))
  expect_error(est_changepoints(data, "ok", nchangepoints,   
               timename, weights, control))
  expect_error(est_changepoints("ok", formula, nchangepoints,   
               timename, weights, control))
})

test_that("check est_regressors", {
  set.seed(1)
  rhos <- est_changepoints(data, formula, nchangepoints,   
               timename, weights, control)
  etas <- est_regressors(rhos, data, formula, timename, weights, control)
  set.seed(1)
  rhos2 <- est_changepoints(data, formula, nchangepoints = 2,  
                timename, weights, list(nit = 20, seed = 1))
  etas2 <- est_regressors(rhos2, data, formula, timename, weights,
                   list(nit = 20, seed = 1))

  expect_is(etas, "matrix")
  expect_equal(colnames(etas), c("1_2:(Intercept)", "2_2:(Intercept)"))
  expect_equal(dim(etas), c(20, 2))
  expect_equal(round(sum(etas[ , 1]), 1), 35.7)
  expect_equal(round(sum(etas2[1:10 , 1]), 1), 15.2)
  expect_error(est_regressors("ok", data, formula, timename,weights,
                control))
  expect_error(est_regressors(rhos, data, formula, weights,timename, "ok"))
  expect_error(est_regressors(rhos, data, formula, "ok",  weights, control))
  expect_error(est_regressors(rhos, data, timename, "ok", 
               timename, control))
  expect_error(est_regressors(rhos, "ok", formula, timename, weights,control))
  expect_error(est_regressors(rhos, data, "ok", timename, weights,control))
  rhosx <- rhos
  names(rhosx)[1] <- "ok"
  expect_error(est_regressors(rhosx, data, formula, timename, weights,
               control))
})

test_that("check package_TS", {
  summ <- package_TS(data, formula, timename, weights,  control, 
                       rho_dist, eta_dist)
  expect_is(summ, "TS_fit")
  expect_equal(length(summ), 17)
  expect_equal(names(summ)[3], "nchangepoints")
  expect_error(
    package_TS("ok", formula, timename, weights, control, rho_dist, 
                  eta_dist))
  expect_error(package_TS(data, formula, "ok", weights,control, rho_dist,
                            eta_dist))
  expect_error(package_TS(data, "ok", timename, weights, control, rho_dist,
                            eta_dist))

  expect_error(package_TS(data, formula, timename, weights,"ok", rho_dist,
                            eta_dist))
  expect_error(package_TS(data, formula, timename,weights, control, "ok",
                            eta_dist))
  expect_error(package_TS(data, formula, timename,weights, control,
                            rho_dist, "ok"))
  expect_error(package_TS("ok", formula, timename,weights, control,
                            rho_dist, eta_dist))
})

test_that("check TS", {
  expect_is(TSmod, "TS_fit")
  expect_equal(length(TSmod), 17)
  expect_equal(TSmod$nchangepoints, 1)
  expect_error(TS(data, formula, nchangepoints = 0, timename, weights, "ok"))
  expect_error(TS(data, formula, nchangepoints = 0, "ok", weights, control))
  expect_error(TS(data, "ok", nchangepoints = 0, timename, "ok", control))
  expect_error(TS("ok", formula, nchangepoints = 0, timename, weights, 
                   control))
  expect_error(TS(data, formula, "ok", timename, weights, control))
})

test_that("check check_TS_inputs", {
  expect_silent(check_TS_inputs(data, formula, nchangepoints, 
                                timename, weights, control))
  expect_error(check_TS_inputs(data, formula, nchangepoints,  
                               timename, weights, "ok"))
  expect_error(check_TS_inputs(data, formula, nchangepoints,
                               timename, "ok", control))
  expect_error(check_TS_inputs(data, formula, "ok",   
                               timename, weights, control))
  expect_error(check_TS_inputs(data, "ok", nchangepoints,   
                               timename, weights, control))
  expect_error(check_TS_inputs("ok", formula, nchangepoints,   
                               timename, weights, control))
  expect_error(check_TS_inputs(data, formula, nchangepoints,   
                               "ok",weights,  control))
})

test_that("check print for TS_fit", {
  expect_output(print(TSmod))
})

test_that("check progress bar functions", {
  expect_is(prep_pbar(), "progress_bar")
  expect_silent(prep_pbar(list(quiet = TRUE)))
  expect_error(prep_pbar("ok"))
  expect_error(prep_pbar(list(), "ok"))
  expect_error(prep_pbar(nr = 0.5))
  pp <- prep_pbar()
  expect_is(update_pbar(pp), "progress_bar")
  expect_silent(update_pbar(pp, list(quiet = TRUE)))

  expect_error(update_pbar(pp, "ok"))
  expect_error(update_pbar("ok", list()))
})


test_that("check AIC for TS_fit", {
  expect_equal(round(AIC(TSmod), 1), 404.2)
})

test_that("check check_formula", {
  expect_silent(check_formula(data, formula))
  expect_error(check_formula(data, ~1))
  expect_error(check_formula(data, ~ok))
  expect_error(check_formula(document_covariate_table, formula))
})

test_that("check TS_control", {
  expect_is(TS_control(), "list")
  expect_equal(length(TS_control()), 16)
})
