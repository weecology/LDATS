
#' @title Simulate LDA data from an LDA structure given parameters
#'
#' @description For a given set of parameters \code{alpha} and \code{Beta} and
#'   document-specific total word counts, simulate a document-by-term matrix.
#'   Additional structuring variables (the numbers of topics (k),
#'   documents (M), terms (V)) are inferred from input objects.
#'
#' @param N A vector of document sizes (total word counts). Must be integer
#'   conformable. Is used to infer the total number of documents.
#' 
#' @param Beta \code{matrix} of categorical distribution parameters defining
#'   terms within topics. Dimension: k x V (number of topics x number of 
#'   terms). Used to infer both (k) and (V). Must be non-negative and sum to 
#'   1 within topics.
#' 
#' @param alpha Single positive numeric value for the Dirichlet distribution
#'   parameter defining topics within documents. To specifically define
#'   document topic probabilities, use \code{Theta}.
#' 
#' @param Theta \code{matrix} of probabilities defining topics within 
#'   documents. Dimension: M x k (documents x topics). Must be non-negative 
#'   and sum to 1 within documents. To generally define document topic 
#'   probabilities, use \code{alpha}.
#' 
#' @param seed Input to \code{\link{set.seed}}.
#'
#' @return A document-by-term \code{matrix} of counts (dim: M x V).
#' 
#' @examples
#'   N <- c(10, 22, 15, 31)
#'   alpha <- 1.2
#'   Beta <- matrix(c(0.1, 0.1, 0.8, 0.2, 0.6, 0.2), 2, 3, byrow = TRUE)
#'   sim_LDA_data(N, Beta, alpha = alpha)
#'   Theta <- matrix(c(), 4, 2
#'   sim_LDA_data(N, Beta, Theta = Theta)
#'
#' @export
#'
sim_LDA_data <- function(N, Beta, alpha = NULL, Theta = NULL, seed = NULL){
  if (length(dim(N)) > 1 | !is.numeric(N) || !all(N %% 1 == 0)){
    stop("N must be a vector of integer conformable values")
  }
  if (is.null(alpha) & is.null(Theta)){
    stop("alpha or Theta must be defined")
  } 
  if (!is.null(alpha) &&
      (length(alpha) != 1 | !is.numeric(alpha) || alpha <= 0)){
    stop("alpha must be a single postitive number")
  }
  if (!is.null(Theta) &&
     (length(dim(Theta)) != 2 | !is.numeric(Theta))){
    stop("Theta must be a 2-dimensional numeric matrix")
  }
  if (!is.null(Theta) &&
     (any(round(rowSums(Theta), 10) != 1) | any(Theta < 0))){
    stop("Theta values must be non-negative and sum to 1 within documents")
  }
  if (length(dim(Beta)) != 2 | !is.numeric(Beta)){
    stop("Beta must be a 2-dimensional numeric matrix")
  }
  if (any(round(rowSums(Beta), 10) != 1) | any(Beta < 0)){
    stop("Beta values must be non-negative and sum to 1 within documents")
  }
  set.seed(seed)
  M <- length(N)
  k <- dim(Beta)[1]
  V <- dim(Beta)[2]
  if(is.null(Theta)){
    if(k == 1){
      Theta <- matrix(1, nrow = M, ncol = 1)
    } else{
      alpha2 <- rep(alpha, k)
      Theta <- rdirichlet(M, alpha2)
    }
  }
  z <- matrix(NA, nrow = M, ncol = k)
  w <- matrix(NA, nrow = M, ncol = V)
  for(d in 1:M){
    z[d, ] <- table(rcat(N[d], Theta[d, ], 1:k))
    wk <- matrix(NA, nrow = k, ncol = V)
    for(ik in 1:k){
      wk[ik,] <- table(rcat(z[d, ik], Beta[ik,], 1:V))
    }
    w[d,] <- apply(wk, 2, sum)
  }
  w
}

#' @title Simulate TS data from a TS model structure given parameters
#'
#' @description For a given set of covariates \code{X}; parameters \code{Eta},
#'   \code{rho}, and \code{err}; and document-specific time stamps \code{tD},
#'   simulate a document-by-topic matrix. Additional structuring variables 
#'   (numbers of topics (k), documents (M), segments (S), and 
#'   covariates per segment (C)) are inferred from input objects.
#'
#' @param X \code{matrix} of covariates, dimension M (number of documents) x 
#'   SC (number of segments x number of covariates, including the intercept).
#' 
#' @param Eta \code{matrix} of regression parameters across the segments,
#'   dimension: SC (number of segments x number of covariates, including the
#'   intercept) x k (number of topics).
#' 
#' @param rho Vector of integer-conformable time locations of changepoints or 
#'   \code{NULL} if no changepoints. Used to determine the number of 
#'   segements. Must exist within the bounds of the times of the documents,
#'   \code{tD}.
#'
#' @param tD Vector of integer-conformable times of the documents. Must be
#'   of length M (as determined by \code{X}). 
#' 
#' @param err Additive error on the link-scale. Must be a non-negative 
#'   \code{numeric} value. Default value of \code{0} indicates no error.
#' 
#' @param seed Input to \code{\link{set.seed}}.
#'
#' @return A document-by-topic \code{matrix} of probabilities (dim: M x k).
#' 
#' @examples
#'   tD <- c(1, 3, 4, 6)
#'   rho <- 3
#'   X <- matrix(c(1,1,0,0,1,2,0,0,0,0,1,3,0,0,1,4), 
#'          nrow = length(tD), ncol = 4, byrow = TRUE)
#'   Eta <- matrix(c(0.5, 1.2, 0.3, 1.1, 0.9, 0.1, 0.5, 0.5), 
#'           nrow = ncol(X), ncol = 2, byrow = TRUE)
#'   sim_TS_data(X, Eta, rho, tD, err = 1)
#'   
#' @export
#'
sim_TS_data <- function(X, Eta, rho, tD, err = 0, seed = NULL){

  if (length(dim(tD)) > 1 | !is.numeric(tD) || !all(tD %% 1 == 0)){
    stop("tD must be a vector of integer conformable values")
  }
  if (!is.null(rho) &&
      (length(dim(rho)) > 1 | !is.numeric(rho) || !all(rho %% 1 == 0))){
    stop("rho must be NULL or a vector of integer conformable values")
  }
  if (length(dim(X)) != 2 | !is.numeric(X)){
    stop("X must be a 2-dimensional numeric matrix")
  }
  if (length(dim(Eta)) != 2 | !is.numeric(Eta)){
    stop("Eta must be a 2-dimensional numeric matrix")
  }
  if (length(err) != 1 | !is.numeric(err) || err  < 0){
    stop("err must be a single non-negative number")
  }
  set.seed(seed)
  P <- length(rho)
  S <- 1 + P
  C <- nrow(Eta) / S
  s_start <- c(min(tD), rho + 1)
  s_end <- c(rho, max(tD))
  EGamma <- matrix(NA, nrow = nrow(X), ncol = ncol(Eta))
  for(s in 1:S){
    in1 <- which(tD >= s_start[s] & tD <= s_end[s])
    in2 <- (C * (s - 1) + 1):(C * S)
    X_Eta <- X[in1,in2] %*% Eta[in2,]
    eps <- rnorm(length(X_Eta), 0, err)
    EGamma[in1,] <- softmax(X_Eta + eps) 
  }
  EGamma
}


#' @title Simulate LDA_TS data from LDA and TS model structures and parameters
#'
#' @description For a given set of covariates \code{X}; parameters 
#'   \code{Beta}, \code{Eta}, \code{rho}, and \code{err}; and 
#'   document-specific time stamps \code{tD} and lengths \code{N}),
#'   simulate a document-by-topic matrix.
#'   Additional structuring variables (the numbers of topics (k), terms (V), 
#'   documents (M), segments (S), and covariates per segment (C))
#'   are inferred from input objects.
#'
#' @param N A vector of document sizes (total word counts). Must be integer
#'   conformable. Is used to infer the total number of documents.
#' 
#' @param Beta \code{matrix} of categorical distribution parameters defining
#'   terms within topics. Dimension: k x V (number of topics x number of 
#'   terms). Used to infer both (k) and (V). Must be non-negative and sum to 
#'   1 within topics.
#'
#' @param X \code{matrix} of covariates, dimension M (number of documents) x 
#'   SC (number of segments x number of covariates, including the intercept).
#' 
#' @param Eta \code{matrix} of regression parameters across the segments,
#'   dimension: SC (number of segments x number of covariates, including the
#'   intercept) x k (number of topics).
#' 
#' @param rho Vector of integer-conformable time locations of changepoints or 
#'   \code{NULL} if no changepoints. Used to determine the number of 
#'   segements. Must exist within the bounds of the times of the documents,
#'   \code{tD}.
#'
#' @param tD Vector of integer-conformable times of the documents. Must be
#'   of length M (as determined by \code{X}). 
#' 
#' @param err Additive error on the link-scale. Must be a non-negative 
#'   \code{numeric} value. Default value of \code{0} indicates no error.
#' 
#' @param seed Input to \code{\link{set.seed}}.
#'
#' @return A document-by-term \code{matrix} of counts (dim: M x V).
#' 
#' @examples
#'   N <- c(10, 22, 15, 31)
#'   tD <- c(1, 3, 4, 6)
#'   rho <- 3
#'   X <- matrix(c(1,1,0,0,1,2,0,0,0,0,1,3,0,0,1,4), 
#'            nrow = length(tD), ncol = 4, byrow = TRUE)
#'   Eta <- matrix(c(0.5, 1.2, 0.3, 1.1, 0.9, 0.1, 0.5, 0.5), 
#'            nrow = ncol(X), ncol = 2, byrow = TRUE)
#'   Beta <- matrix(c(0.1, 0.1, 0.8, 0.2, 0.6, 0.2), 2, 3, byrow = TRUE)
#'   err <- 1
#'   sim_LDA_TS_data(N, Beta, X, Eta, rho, tD, err)
#'   
#' @export
#'
sim_LDA_TS_data <- function(N, Beta, X, Eta, rho, tD, err = 0, seed = NULL){
  EGamma <- sim_TS_data(X, Eta, rho, tD, err, seed)
  sim_LDA_data(N, Beta, Theta = EGamma, seed = seed)
}


