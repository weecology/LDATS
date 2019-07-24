context("Check ptMCMC functions")

# use old RNG method for sample (for test reproducibility)
if ("sample.kind" %in% names(formals(RNGkind))){
  suppressWarnings(RNGkind(sample.kind = "Rounding"))
}

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
topics <- 2
nseeds <- 1
formulas <- ~ 1
nchangepoints <- 1
weights <- document_weights(document_term_table)
timename <- "newmoon"
LDAs <- LDA_set(document_term_table, topics, nseeds)
LDA_models <- select_LDA(LDAs)
control <- list(nit = 20, seed = 1)
mods <- expand_TS(LDA_models, formulas, nchangepoints)
formula <- mods$formula[[1]]
nchangepoints <- mods$nchangepoints[1]
data <- prep_TS_data(document_covariate_table, LDA_models, mods, 1)

set.seed(1)
rho_dist0 <- est_changepoints(data, formula, nchangepoints = 0, timename, 
                              weights, control)
rho_dist <- est_changepoints(data, formula, nchangepoints, timename, 
                              weights, control)
eta_dist <- est_regressors(rho_dist, data, formula, timename, weights, 
                           control)


saves <- prep_saves(nchangepoints, control)
inputs <- prep_ptMCMC_inputs(data, formula, nchangepoints, timename, weights, 
                             control)
cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
ids <- prep_ids(control)

test_that("check prep_proposal_dist", {
  pd <- prep_proposal_dist(nchangepoints, control)
  expect_equal(length(pd), 2)
  expect_equal(dim(pd[[1]]), c(20, 6))
  pd2 <- prep_proposal_dist(0, control)
  expect_equal(length(pd2), 2)
  expect_equal(dim(pd2[[1]]), c(20, 6))
  expect_error(prep_proposal_dist("ok", control))
  expect_error(prep_proposal_dist(nchangepoints, "ok"))
})



test_that("check prep_ptMCMC_inputs", {
  inpts <- prep_ptMCMC_inputs(data, formula, nchangepoints, timename, weights, 
                              control)
  expect_is(inpts, "ptMCMC_inputs")
  expect_equal(length(inpts[[2]]), 6)

  expect_error(
    prep_ptMCMC_inputs("ok", formula, nchangepoints, timename, weights,  
                        control))
  expect_error(
    prep_ptMCMC_inputs(data, "ok", nchangepoints, timename, weights, control))
  expect_error(
    prep_ptMCMC_inputs(data, formula, "ok", timename, weights, control))
  expect_error(
    prep_ptMCMC_inputs(data, formula, nchangepoints, timename, "ok", control))
  expect_error(
    prep_ptMCMC_inputs(data, formula, nchangepoints, "ok", weights, control))
  expect_error(suppressWarnings(
    prep_ptMCMC_inputs(data, formula, nchangepoints, timename, weights,  
                       "ok")))  
})



test_that("check prep_ids", {
  expect_equal(prep_ids(TS_control()), 1:6)
  expect_error(prep_ids("ok"))
  expect_error(prep_ids(list(ntemps = 0.3)))
})
test_that("check update_ids", {
  set.seed(123)
  steps <- step_chains(1, cpts, inputs)
  swaps <- swap_chains(steps, inputs, ids)
  ids <- update_ids(ids, swaps)
  expect_equal(ids, c(1, 2, 4, 3, 5, 6))
})


test_that("check proposed_step_mods", {
  pdist <- inputs$pdist
  ntemps <- length(inputs$temps)
  selection <- cbind(pdist$which_steps[1, ], 1:ntemps)
  prop_changepts <- cpts$changepts
  curr_changepts_s <- cpts$changepts[selection]
  prop_changepts_s <- curr_changepts_s + pdist$steps[1, ]
  if(all(is.na(prop_changepts_s))){
    prop_changepts_s <- NULL
  }
  prop_changepts[selection] <- prop_changepts_s
  mods <- proposed_step_mods(prop_changepts, inputs)

  expect_is(mods, "list")
  expect_is(mods[[1]], "multinom_TS_fit")
  expect_is(mods[[1]][[1]], "list")
  expect_is(mods[[1]][[1]][[1]], "multinom")
  expect_equal(round(mods[[1]][[1]][[1]]$deviance, 1), 43.9)
})

test_that("check propose_step", {
  prop_step <- propose_step(1, cpts, inputs)
  expect_equal(length(prop_step), 2)
  expect_equal(names(prop_step), c("changepts", "lls"))
  expect_equal(prop_step[[1]][1,1], 198)
})
test_that("check eval_step", {
  set.seed(1)
  prop_step <- propose_step(1, cpts, inputs)
  accept_step <- eval_step(1, cpts, prop_step, inputs)
  expect_equal(accept_step, c(T, F, T, T, T, T))
})
test_that("check take_step", {
  prop_step <- propose_step(1, cpts, inputs)
  accept_step <- eval_step(1, cpts, prop_step, inputs)
  taken <- take_step(cpts, prop_step, accept_step)
  expect_equal(length(taken), 3)
  expect_equal(names(taken), c("changepts", "lls", "accept_step"))
  expect_equal(taken[[3]][3], TRUE)
})


test_that("check step_chains", {
  steps <- step_chains(1, cpts, inputs)
  expect_equal(length(steps), 3)
  expect_equal(names(steps), c("changepts", "lls", "accept_step"))
  expect_equal(steps[[3]][3], TRUE)
})

test_that("check swap_chains", {
  steps <- step_chains(1, cpts, inputs)
  swaps <- swap_chains(steps, inputs, ids)
  expect_equal(length(swaps), 4)
  expect_equal(names(swaps), c("changepts", "lls", "ids", "accept_swap"))
  expect_equal(swaps[[3]][3], 4)
})

test_that("check count_trips", {
  set.seed(1)
  expect_equal(length(count_trips(rho_dist$ids)), 2)
  expect_equal(names(count_trips(rho_dist$ids)), 
               c("trip_counts", "trip_rates"))
  expect_equal(count_trips(rho_dist$ids)[[1]][[3]], 0)

  expect_equal(count_trips(matrix(1, 6, 100))[[1]][3], 0)
})

test_that("check diagnose_ptMCMC", {
  set.seed(1)
  expect_equal(diagnose_ptMCMC(rho_dist0), NULL)
  expect_equal(length(diagnose_ptMCMC(rho_dist)), 4)
  expect_equal(names(diagnose_ptMCMC(rho_dist)), 
    c("step_acceptance_rate", "swap_acceptance_rate", "trip_counts",
      "trip_rates"))
  expect_equal(diagnose_ptMCMC(rho_dist)[[1]][1], 0.25)
})

test_that("check prep_saves", {
  svs <- prep_saves(nchangepoints, control)
  expect_is(svs, "list")
  expect_equal(length(svs), 5)
  expect_equal(dim(svs[[1]]), c(1, 6, 20))
  expect_error(prep_saves("ok", control))
  expect_error(prep_saves(nchangepoints, "ok"))
})
test_that("check update_saves", {
  set.seed(1)
  saves <- prep_saves(nchangepoints, control)
  inputs <- prep_ptMCMC_inputs(data, formula, nchangepoints,  
                               timename, weights, control)
  cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
  ids <- prep_ids(control)
  steps <- step_chains(1, cpts, inputs)
  swaps <- swap_chains(steps, inputs, ids)
  saves <- update_saves(1, saves, steps, swaps)
  expect_is(saves, "list")
  expect_equal(length(saves), 5)
  expect_equal(dim(saves[[1]]), c(1, 6, 20))
  expect_equal(saves[[1]][1, 1, 1], 309)
})

test_that("check process_saves", {
  set.seed(1)
  saves <- prep_saves(nchangepoints, control)
  inputs <- prep_ptMCMC_inputs(data, formula, nchangepoints,  
                               timename, weights, control)
  cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
  ids <- prep_ids(control)
  for(i in 1:control$nit){
    steps <- step_chains(i, cpts, inputs)
    swaps <- swap_chains(steps, inputs, ids)
    saves <- update_saves(i, saves, steps, swaps)
    cpts <- update_cpts(cpts, swaps)
    ids <- update_ids(ids, swaps)
  }
  out <- process_saves(saves, control)
  expect_is(out, "list")
  expect_equal(length(out), 5)
  expect_equal(dim(out[[1]]), c(1, 6, 20))
  expect_equal(out[[1]][1, 1, 1], 309)
  expect_equal(out[[1]][1, 1, 20], 270)
  out2 <- process_saves(saves, list(burnin = 10, nit = 20))
  expect_is(out2, "list")
  expect_equal(length(out2), 5)
  expect_equal(dim(out2[[1]]), c(1, 6, 10))
  expect_equal(out2[[1]][1, 1, 1], 309)
  expect_equal(out2[[1]][1, 1, 3], 301)
})

test_that("check prep_cpts", {
  set.seed(1)
  cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
  expect_is(cpts, "list")
  expect_equal(length(cpts), 2)
  expect_equal(cpts[[1]][1,1], 268)

  expect_error(prep_cpts("ok", formula, nchangepoints, timename, weights, 
                         control))
  expect_error(prep_cpts(data, "ok", nchangepoints, timename, weights, 
                         control))
  expect_error(prep_cpts(data, formula, "ok", timename, weights, control))
  expect_error(prep_cpts(data, formula, nchangepoints, "ok", weights, 
                         control))
  expect_error(prep_cpts(data, formula, nchangepoints, timename, "ok", 
                         control))
  expect_error(prep_cpts(data, formula, nchangepoints, timename, weights,  
                         "ok"))
})
test_that("check update_cpts", {
  set.seed(1)
  saves <- prep_saves(nchangepoints, control)
  inputs <- prep_ptMCMC_inputs(data, formula, nchangepoints, timename, 
                               weights, control)
  cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
  ids <- prep_ids(control)
  steps <- step_chains(1, cpts, inputs)
  swaps <- swap_chains(steps, inputs, ids)
  cpts <- update_cpts(cpts, swaps)
  expect_is(cpts, "list")
  expect_equal(length(cpts), 2)
  expect_equal(cpts[[1]][1,1], 309)  
})

test_that("check prep_temp_sequence", {
  expect_equal(length(prep_temp_sequence()), 6)
  expect_equal(length(prep_temp_sequence(list(ntemps = 9))), 9)
  expect_equal(round(prep_temp_sequence()[3], 2), 8)
  expect_equal(round(prep_temp_sequence(list(q = 1))[3], 1), 2.8)
  expect_error(prep_temp_sequence(123))
})

