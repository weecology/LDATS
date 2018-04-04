  chunk_memo <- memoise::memoise(multinom_chunk)


system.time(
for(W in 1:1000){
xx<- chunk_memo(formula_RHS, data, start_time = start_times[i], 
                   end_time = end_times[i], weights)

}
)

system.time(
for(W in 1:1000){
xx<- multinom_chunk(formula_RHS, data, start_time = start_times[i], 
                   end_time = end_times[i], weights)

}
)



system.time(
for(W in 1:1000){
xx<- ts_memo(formula_RHS, data, changepoints , weights)
}
)