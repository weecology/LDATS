
TS <- function(gamma, formula, nchangepoints, weights, 
               ptMCMC_controls = ptMCMC_controls_list()){



}

ptMCMC_controls_list <- function(){
  out <- list()
  class(out) <- c("ptMCMC_controls", "list")
  out
}