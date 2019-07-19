devtools::load_all()

# exploring here just some basic simulation and recovery with the LDA

  set.seed(123)
  N <- sample(100:100, 1000, TRUE)
  M <- length(N)  
  k <- 3
  V <- 200
  Beta <- matrix(NA, k, V)
  for(i in 1:k){
    Beta[i,] <- rep(1/V, V)
  }

  alpha <- 1.2
  Theta <- matrix(NA, M, k)
  for(d in 1:M){
    Theta[d,] <- rep(1/k, k)
  }

  simmed_a <- sim_LDA_data(N, Beta, alpha = alpha)
  simmed_t <- sim_LDA_data(N, Beta, Theta = Theta)

  fit_a <- LDA(simmed_a, 3)
  fit_a@alpha
  fit_t <- LDA(simmed_t, 3)
  fit_t@gamma

  set_a <- LDA_set(simmed_a, 2:5, 100)
  select_LDA(set_a)
  set_t <- LDA_set(simmed_t, 2:5, 100)
  select_LDA(set_t)

