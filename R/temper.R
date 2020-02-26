

# use the existing temperature schematic as a start
#  build a template into the temp_fun

# now need to generalize prep_temp_sequence, but leverages existing control
#  list approach exactly how we want


temp_fun <- function(chain = 1, control = list()){
  prep_temp_sequence(control = control)[chain]
}



