# This script assumes that `dat` is defined

library(topicmodels)


#' Compute hellinger distance
#' 
#' Hellinger distance is Euclidean distance between square roots, divided by
#' sqrt(2)
#' 
#' @param a,b numeric vectors
#' 
Hellinger = function(a, b){
  diff = sqrt(a) - sqrt(b)
  sqrt(sum(diff^2) / 2)
}

#' Find the configuration with the minimum Hellinger distance between two
#' probability matrices like those from `ps`.
#' 
#' @param p1,p2 numeric vectors: species composition
#' 
min_H = function(p1, p2) {
  # Find the cost associated with each pairwise topic assignment
  costs = matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      costs[i, j] = Hellinger(p1[i, ], p2[j, ])
    }
  }
  # Minimize the sum of the pairwise costs.
  # Adding a constant doesn't change the results, but all costs must be strictly
  # positive.
  assignment = clue::solve_LSAP(1 + costs)
  list(
    min_cost = mean(costs[cbind(1:k, assignment)]),
    assignment = assignment
  )
}


# seed used to make the figures in rodent_LDA_analysis.R
SEED = 113052032
 #SEED = 28 has a lower loglik

#'
#'
#'
#' @param ldas list of results from LDA model
#' @param seedlist vector of seeds for which LDA models will be fit and results compared
#'
calculate_LDA_distance = function(ldas) {
  
  # Log-likelihoods for each lda
  lls = purrr:::map_dbl(ldas, logLik)
  
  # Pick the LDA model with the highest log-likelihood, as is done in the main
  # analysis.
  best_seed = seedlist[which.max(lls)]
  
  # Discard models whose log-likelihood is much lower than the "best" model
  #ldas = ldas[lls > max(lls) - 100]
  
  # Topic allocations to each species (probabilities that sum to 1 in each row)
  ps = purrr::map(ldas, ~exp(.x@beta))
  ps_best = exp(ldas[1]$beta)
  
  # Calculate Hellinger distances from this model and find the farthest one
  minimum_distances = lapply(
    1:length(ps),
    function(i){min_H(ps[[best_lda]], ps[[i]])})
  farthest_lda = which.max(purrr::map_dbl(minimum_distances, "min_cost"))
  
  # Find the species allocations of each topic for these two models
  best_ps = ps[[best_lda]]
  farthest_ps = ps[[farthest_lda]][minimum_distances[[farthest_lda]]$assignment, ]
  
  # Re-order the rows of the "farthest" model so it matches the "best" model
  assignment = min_H(best_ps, farthest_ps)$assignment
  

}




#' plot results
#' "best" is in black, "farthest" is in red
#' 
plot_spcomp_distance = function(dat,best_ps,farthest_ps) {
  par(mfrow = c(2, 2))
  for (row in 1:k) {
    plot(NULL, xlim = c(1, 21), ylim = c(0, 1),
         axes = FALSE, main = paste("topic", row), ylab = "", xlab = "",
         yaxs = "i")
    abline(v = 1:ncol(dat), col = scales::alpha(1, .05))
    points(best_ps[row, ], col = 1, pch = 16)
    points(farthest_ps[assignment[row], ], col = 2, lwd = 2)
    axis(2)
    axis(1, 1:ncol(dat), colnames(dat), las = 2)
  }
  par(mfrow = c(1, 1))
  
}

