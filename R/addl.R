
#
# using this script to help tidy up the time_series_functions script
#  functions here still need documentation etc etc

prep_formula <- function(formula){
  character_formula <- as.character(formula)
  character_formula[length(character_formula)]
}

memoise_TF <- function(function_name = "multinom_ts", memo = TRUE){
  function_entity <- eval(parse(text = function_name))
  if(memo){
    out <- memoise(function_entity)
  } else{
    out <- function_entity
  }
  out
}