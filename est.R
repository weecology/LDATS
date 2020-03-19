

est_changepointsx <- function(data, formula, nchangepoints, timename, weights, 
                             control = list()){
  check_TS_inputs(data, formula, nchangepoints, timename, weights, control)
  control <- do.call("TS_control", control)
  data <- time_order_data(data, timename = timename)
  if (nchangepoints == 0){
    return(NULL)
  }

# break this into a "classic" approach or whatever and give it its own
#  function that is akin to temper

  saves <- prep_saves(nchangepoints, control)

  inputs <- prep_ptMCMC_inputs(data, formula, nchangepoints, timename, 
                               weights, control)
  cpts <- prep_cpts(data, formula, nchangepoints, timename, weights, control)
  ids <- prep_ids(control)
  pbar <- prep_pbar(control, "rho")

  for(i in 1:control$nit){
    update_pbar(pbar, control)
    steps <- step_chains(i, cpts, inputs)
    swaps <- swap_chains(steps, inputs, ids)
    saves <- update_saves(i, saves, steps, swaps)
    cpts <- update_cpts(cpts, swaps)
    ids <- update_ids(ids, swaps)
  }
  process_saves(saves, control)

# the separate function runs to here
# and then process the output


}
