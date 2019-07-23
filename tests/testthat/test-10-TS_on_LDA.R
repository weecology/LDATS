context("Check TS_on_LDA functions")

data(rodents)
lda_data <- rodents$document_term_table
document_term_table <- rodents$document_term_table
document_covariate_table <- rodents$document_covariate_table
topics <- 2
nseeds <- 1
formulas <- ~ 1
nchangepoints <- 0
weights <- document_weights(document_term_table)
timename <- "newmoon"
LDAs <- LDA_set(document_term_table, topics, nseeds)
LDA_models <- select_LDA(LDAs)


test_that("check prep_TS_data", {
  mods <- expand_TS(LDA_models, formulas=c(~1, ~sin_year, ~newmoon), 
                    nchangepoints = 0:1)
  nmods <- nrow(mods)
  TSmods <- vector("list", nmods)

  for (i in 1:nmods){
    data_i <- prep_TS_data(document_covariate_table, LDA_models, mods, i)
    expect_is(data_i, "data.frame")
    expect_equal(dim(data_i), c(436, 4))
  }

  expect_error(prep_TS_data(document_covariate_table, LDA_models, mods, "ok"))
  expect_error(prep_TS_data(document_covariate_table, LDA_models, "ok", i))
  expect_error(prep_TS_data(document_covariate_table, "ok", mods, i))
  expect_error(prep_TS_data("ok", LDA_models, mods, i))

  data_1 <- prep_TS_data(document_covariate_table, LDA_models[[1]], mods, 1)
  data_i <- prep_TS_data(document_covariate_table, LDA_models, mods, 1)
  expect_equal(data_1, data_i)
})

test_that("check select_TS", {
  mods <- TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  nchangepoints = 0:1, timename, weights, 
                  control = list(nit = 1e2))
  sel_mod <- select_TS(mods, list())

  expect_equal(length(mods), 2)
  expect_equal(length(sel_mod), 17)
  expect_is(mods, "TS_on_LDA")
  expect_is(sel_mod, "TS_fit")
  expect_error(select_TS(mods, "ok"))
  expect_error(select_TS("ok", list()))

  xxfun <- function(x){x}
  expect_warning(select_TS(mods, list(selector = xxfun)))
})

test_that("check package_TS_on_LDA", {
  mods <- expand_TS(LDA_models, formulas, nchangepoints = 0:1)
  nmods <- nrow(mods)
  TSmods <- vector("list", nmods)

  for(i in 1:nmods){
    print_model_run_message(mods, i, LDA_models, list())
    formula_i <- mods$formula[[i]]
    nchangepoints_i <- mods$nchangepoints[i]
    data_i <- prep_TS_data(document_covariate_table, LDA_models, mods, i)
    TSmods[[i]] <- TS(data_i, formula_i, nchangepoints_i,timename, weights, 
                      control = list(nit = 1e2))
  }
  expect_is(package_TS_on_LDA(TSmods, LDA_models, mods), "TS_on_LDA")
  expect_is(package_TS_on_LDA(TSmods, LDA_models[[1]], mods), "TS_on_LDA")
  expect_error(package_TS_on_LDA(TSmods, LDA_models, "ok"))
  expect_error(package_TS_on_LDA("ok", LDA_models, mods))
  expect_error(package_TS_on_LDA(TSmods, "ok", mods))
})


test_that("check TS_on_LDA", {
  mods <- TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  nchangepoints = 0:1, timename, weights, 
                  control = list(nit = 1e2))
  expect_is(mods, "TS_on_LDA")
  expect_equal(length(mods), 2)
  expect_is(mods[[1]], "TS_fit")
  expect_is(mods[[2]], "TS_fit")

  expect_error(TS_on_LDA())
  expect_error(TS_on_LDA("ok", document_covariate_table, formulas, 
                  nchangepoints, timename, weights, list()))
  expect_error(TS_on_LDA(LDA_models, "ok", formulas, 
                  nchangepoints, timename, weights, list()))
  expect_error(TS_on_LDA(LDA_models, document_covariate_table, "ok", 
                  nchangepoints, timename, weights, list()))
  expect_error(TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  "ok", timename, weights, list()))
  expect_error(TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  nchangepoints, timename, "ok", list()))
  expect_error(TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  nchangepoints, "ok", weights, list()))
  expect_error(TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  nchangepoints, timename, weights, "ok"))
})


test_that("check printing for TS_on_LDA", {
  m1 <- TS_on_LDA(LDA_models, document_covariate_table, formulas, 
                  nchangepoints, timename,weights, 
                  control = list(nit = 1e2))
  expect_output(print(m1))
})

test_that("check print_model_run_message", {
  mods <- expand_TS(LDA_models, formulas, nchangepoints)
  expect_message(print_model_run_message(mods, 1, LDA_models, list()))
  expect_silent(print_model_run_message(mods, 1, LDA_models, 
                          control = list(quiet = TRUE)))
})

test_that("check expand_TS", {
  exp_TS <- expand_TS(LDA_models, formulas, nchangepoints)
  expect_is(exp_TS, "data.frame")
  expect_equal(dim(exp_TS), c(1, 3))
  exp_TS <- expand_TS(LDA_models, c(~1, ~newmoon), 3:10)
  expect_is(exp_TS, "data.frame")
  expect_equal(dim(exp_TS), c(16, 3))
  exp_TS <- expand_TS(LDA_models[[1]], c(~1, ~newmoon), 3:10)
  expect_is(exp_TS, "data.frame")
  expect_equal(dim(exp_TS), c(16, 3))
  expect_error(expand_TS("ok", formulas, nchangepoints))
  expect_error(expand_TS(LDA_models, "ok", nchangepoints))
  expect_error(expand_TS(LDA_models, c("~1", "ok"), nchangepoints))
  expect_error(expand_TS(LDA_models, list("~1", "ok"), nchangepoints))
  expect_error(expand_TS(LDA_models, formulas, 2.5))
})

test_that("check check_nchangepoints", {
  expect_silent(check_nchangepoints(1))
  expect_silent(check_nchangepoints(0))
  expect_silent(check_nchangepoints(1:10))
  expect_error(check_nchangepoints(2.5))
  expect_error(check_nchangepoints(-1))
  expect_error(check_nchangepoints(NULL))
})


test_that("check check_weights", {
  expect_equal(check_weights(TRUE), NULL)
  expect_error(check_weights(FALSE))
  expect_silent(check_weights(weights))
  expect_silent(check_weights(1))
  expect_silent(check_weights(NULL))
  expect_error(check_weights("ok"))
  expect_error(check_weights(-1))
  expect_warning(check_weights(100))
})


test_that("check check_LDA_models", {
  expect_silent(check_LDA_models(LDA_models))
  expect_silent(check_LDA_models(LDA_models[[1]]))
  expect_error(check_LDA_models("ok"))
})


test_that("check check_document_covariate_table", {
  expect_silent(check_document_covariate_table(document_covariate_table))
  expect_silent(check_document_covariate_table(document_covariate_table, 
                                               LDA_models = LDA_models))
  expect_silent(check_document_covariate_table(document_covariate_table, 
                                               LDA_models = LDAs))
  expect_silent(check_document_covariate_table(document_covariate_table, 
                                  document_term_table = document_term_table))
  expect_error(check_document_covariate_table(document_covariate_table,
                                  LDA_models = 1))
  expect_error(check_document_covariate_table(document_covariate_table,
                                  document_term_table = 1))
  expect_error(check_document_covariate_table(document_covariate_table = 1,
                                  LDA_models = LDA_models))
  expect_error(check_document_covariate_table(document_covariate_table = 1,
                                  document_term_table = document_term_table))
  expect_error(check_document_covariate_table(lm(1~1), LDA_models))
})

test_that("check check_timename", {
  expect_silent(check_timename(document_covariate_table, timename))
  expect_error(check_timename("ok", timename))
  expect_error(check_timename(document_covariate_table, "ok"))
  expect_error(check_timename(document_covariate_table, 1))
  expect_error(check_timename(document_covariate_table, 
                              rep(timename, 2)))
  expect_error(check_timename(document_covariate_table, 1))
  expect_error(check_timename(data.frame(letters), "letters"))

  dct2 <- document_covariate_table
  dct2[,timename] <- dct2[,timename] + 0.1
  expect_error(check_timename(dct2, timename))
})

test_that("check check_formulas", {
  expect_silent(check_formulas(formulas, document_covariate_table, list()))
  expect_error(check_formulas("ok", document_covariate_table, list()))
  expect_error(check_formulas(~newmoon, "ok", list()))
  expect_error(check_formulas(c(~1, "ok"), 
               document_covariate_table, list()))
  expect_error(check_formulas(list(~1, "ok"), 
               document_covariate_table, list()))
  expect_error(check_formulas(formulas, document_covariate_table, "ok"))
})


test_that("check check_TS_on_LDA_inputs", {
  expect_silent(
    check_TS_on_LDA_inputs(LDA_models, document_covariate_table, formulas, 
                           nchangepoints, timename, weights, list()))
  expect_error(
    check_TS_on_LDA_inputs(LDA_models, document_covariate_table, formulas, 
                           nchangepoints, timename, weights, "ok"))
  expect_error(
    check_TS_on_LDA_inputs(LDA_models, document_covariate_table, formulas, 
                           nchangepoints, timename, "ok", list()))
  expect_error(
    check_TS_on_LDA_inputs(LDA_models, document_covariate_table, formulas, 
                           "ok", weights, timename, weights, list()))
  expect_error(
    check_TS_on_LDA_inputs(LDA_models, document_covariate_table, "ok", 
                           nchangepoints, timename, weights, list()))
  expect_error(
    check_TS_on_LDA_inputs(LDA_models, "ok", formulas, 
                           nchangepoints, timename, weights, list()))
  expect_error(
    check_TS_on_LDA_inputs("ok", document_covariate_table, formulas, 
                           nchangepoints, timename,  weights, list()))
  expect_error(
    check_TS_on_LDA_inputs(LDA_models, document_covariate_table, formulas, 
                           nchangepoints, "ok",  weights, list()))

})
