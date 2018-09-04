

MTSn <- function(data, formula = ~1, nchangepoints = 1, weights = NULL, 
                 memo = TRUE
                 # arguments from here down could go into an MCMC options list
                 ){

  formula <- prep_formula(formula)
  ts_memo <- memoise_TF("multinom_ts", memo)
  change_points <- estimate_cps()
  betas <- estimate_betas()
  MTS_output(match.call(), change_points, betas)
}




estimate_cps <- function(){

}


estimate_betas <- function(){


}

MTS_output <- function(main_call, change_points, betas){


}



