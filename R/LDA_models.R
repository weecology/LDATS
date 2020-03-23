
#' @title Latent Dirichlet Allocation Linguistic Decomposition Analysis
#'   as conducted via the topicmodels package
#'
#' @description Fit the standard LDATS LDA model (a true Latent Dirichlet
#'   Allocation) using \code{\link[topicmodels]{LDA}} (Grun and Hornik 2011).
#'   Default methodology is the Variational Expectation Maximization routine 
#'   (VEM) as described by Blei \emph{et al.} (2003) and implemented by 
#'   Grun and Hornik (2011). \cr \cr
#'   If the model is defined to only fit one topic, \code{\link{identity_LDA}}
#'   is used by default.
#'
#' @param LDA A prepared (via \code{\link{prepare_LDA}} LDA model
#'   \code{list}.
#'
#' @param ... Additional arguments to be passed to 
#'   \code{\link[topicmodels]{LDA}} as a \code{control} input.
#'
#' @param seeded \code{logical} indicator of if the LDA should be a seeded
#'   replicate. 
#'
#' @param method Fitting routine used in \code{\link[topicmodels]{LDA}}.
#'   Currenlty, only \code{"VEM"} and \code{"Gibbs"} are supported.
#'
#' @return \code{LDA} \code{list} with components
#'  \describe{
#'    \item{alpha}{parameter estimate.}
#'    \item{beta}{parameter estimate.}
#'    \item{terms}{\code{character} \code{vector} of term names.}
#'    \item{document_topic_matrix}{estimated latent topic compositions.}
#'    \item{test_document_topic_matrox}{estimated latent topic compositions 
#'      of the test data (not presently available for usage).}
#'    \item{log_likelihood}{model log likelihood.}
#'    \item{data}{data object used to fit the LDA model.}
#'    \item{data_subset}{number of the data subset from the whole data set.}
#'    \item{topics}{\code{integer} number of topics in the model.}
#'    \item{replicat}{\code{integer} replicate number.}
#'    \item{control}{\code{list} of controls used to fit the model. See
#'      \code{\link{LDA_control}}.}
#'  }
#' 
#' @references 
#'   Blei, D. M., A. Y. Ng, and M. I. Jordan. 2003. Latent Dirichlet
#'   Allocation. \emph{Journal of Machine Learning Research} 
#'   \strong{3}:993-1022.
#'   \href{http://jmlr.csail.mit.edu/papers/v3/blei03a.html}{link}.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#' 
#' @export
#'
topicmodels_LDA <- function(LDA, method = "VEM", seeded = TRUE, ...){
  data <- LDA$data
  topics <- LDA$topics 
  rep <- LDA$rep
  data_subset <- LDA$data_subset
  if(topics == 1){
    identity_LDA(LDA)
  } else{
    fun_control <- list(...)
    if(seeded){
      fun_control <- update_list(fun_control, seed = rep * 2)
    }
    mod <- topicmodels::LDA(x = data$train$document_term_table, k = topics, 
                            method = method, control = fun_control)
    mod_ll <- sum(mod@loglikelihood)
    alpha <- tryCatch(as.integer(mod@control@estimate.alpha), 
                      error = function(x){0})
    df <- alpha + length(mod@beta)
    attr(mod_ll, "df") <- df
    attr(mod_ll, "nobs") <- mod@Dim[1] * mod@Dim[2]
    class(mod_ll) <- "logLik"
    out <- update_list(LDA, params = list(alpha = mod@alpha, beta = mod@beta),
                            document_topic_table = mod@gamma, 
                            terms = mod@terms,
                            test_document_topic_matrix = NULL, #not available
                            log_likelihood = mod_ll)
    class(out) <- c("LDA", "list")
    out
  } 
}

#' @title Identity Linguistic Decomposition Analysis
#'
#' @description This function acts as an "identity" model, wherein the 
#'   output is functionally the input. This allows for "single-topic" models
#'   that do not actually decompose the data to be included in the model set.
#'
#' @param LDA A prepared (via \code{\link{prepare_LDA}} LDA model
#'   \code{list}.
#'
#' @return \code{LDA} \code{list} with components (many of which are 
#'  placeholders):
#'  \describe{
#'    \item{alpha}{parameter estimate.}
#'    \item{beta}{parameter estimate.}
#'    \item{terms}{\code{character} \code{vector} of term names.}
#'    \item{document_topic_matrix}{estimated latent topic compositions.}
#'    \item{test_document_topic_matrox}{estimated latent topic compositions 
#'      of the test data (not presently available for usage).}
#'    \item{log_likelihood}{model log likelihood.}
#'    \item{data}{data object used to fit the LDA model.}
#'    \item{data_subset}{number of the data subset from the whole data set.}
#'    \item{topics}{\code{integer} number of topics in the model.}
#'    \item{replicat}{\code{integer} replicate number.}
#'    \item{control}{\code{list} of controls used to fit the model. See
#'      \code{\link{LDA_control}}.}
#'  }
#' 
#' @export
#'
identity_LDA <- function(LDA){
  data <- LDA$data
  rep <- LDA$rep
  data_subset <- LDA$data_subset
  document_topic_table <- matrix(1, ncol = 1,
                                 nrow = NROW(data$train$document_term_table))
  beta <- apply(data$train$document_term_table, 2, sum)
  beta <- log(beta/sum(beta))
  beta <- matrix(beta, nrow = 1)
  out <- update_list(LDA,  params = list(alpha = 1, beta = beta), 
                           document_topic_table = document_topic_table, 
                           log_likelihood = NULL, data = data,
                           topics = 1, rep = rep, data_subset = data_subset,
                           terms = colnames(data$train$document_term_table),
                           test_document_topic_matrix = NULL) #not available
  class(out) <- c("LDA", "list")
  out
}
