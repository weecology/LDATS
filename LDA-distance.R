# This script assumes that `dat` is defined

library(topicmodels)

N_seeds = 100
k = 4

# Fit a bunch of LDA models with different seeds
# Only use every other seed because consecutive seeds give identical results (!?)
ldas = purrr::map(2 * seq(N_seeds), 
                  ~LDA(dat, k = 4, method = "VEM", control = list(seed = .x)))

# Topic allocations to each species (probabilities sum to 1)
ps = purrr::map(ldas, ~exp(.x@beta))

# Log-likelihoods
lls = purrr:::map_dbl(ldas, logLik)

# Hellinger distance is Euclidean distance between square roots, divided by
# sqrt(2)
H = function(a, b){
  diff = sqrt(a) - sqrt(b)
  sqrt(sum(diff^2) / 2)
}

# Find the configuration with the minimum Hellinger distance between two
# probability matrices like those from `ps`.
min_H = function(p1, p2) {
  # Find the cost associated with each pairwise topic assignment
  costs = matrix(0, k, k)
  for (i in 1:k) {
    for (j in 1:k) {
      costs[i, j] = H(p1[i, ], p2[j, ])
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


dists = matrix(0, N_seeds, N_seeds)

# Fill in the upper triangle of the distance matrix with the minimum
# distance
for (i in 1:(N_seeds - 1)) {
  for (j in (i+1):N_seeds) {
    result = min_H(ps[[i]], ps[[j]])
    dists[i, j] = result$min_cost
  }
}
dists = dists + t(dists) # Fill in the lower triangle


# We can probably ignore the models with poor likelihoods
include = lls > max(lls) - 50
dists = dists[include, include]

# Pick the LDA model with the highest log-likelihood, as is done in the main
# analysis.
best_lda = which.max(lls[include])

# Pick the LDA model that is farthest from the best model, among our "included"
# models.
farthest_lda = which.max(dists[best_lda, ])

# Find the species allocations of each topic for these two models
best_ps = ps[include][[best_lda]]
farthest_ps = ps[include][[farthest_lda]]

# Re-order the rows of the "farthest" model so it matches the "best" model
assignment = min_H(best_ps, farthest_ps)$assignment



# "best" is in black, "farthest" is in red
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
