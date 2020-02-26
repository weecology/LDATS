
logLik.LDA_VEM <- function(object, ...){
  val <- sum(object@loglikelihood)
  df <- as.integer(object@control@estimate.alpha) + length(object@beta)
  attr(val, "df") <- df
  attr(val, "nobs") <- object@Dim[1] * object@Dim[2]
  class(val) <- "logLik"
  val
}


logLik_topicmodels_LDA <- function(object, method){
  if (method %in% c("VEM", "Gibbs")){
    val <- sum(object@loglikelihood)
    df <- as.integer(object@control@estimate.alpha) + length(object@beta)
    attr(val, "df") <- df
    attr(val, "nobs") <- object@Dim[1] * object@Dim[2]
    class(val) <- "logLik"
    val
  } else if (isnull(method)){
    logLik(object)
  } else{
    stop("method not recognized")
  }
}
