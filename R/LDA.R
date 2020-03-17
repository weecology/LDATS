# LDA is now functionally what LDA_set was before
# "reps" replaces "nseeds" but otherwise things are the same generally at api
# the control list is expanded to now include a base LDA modeling function
# and then arguments for all three functions, it also gets the subset info
# because it may be needed, and in addition to quiet a boolean "soften"
# which softens errors during model running
# prep_LDA_models combines the data (conforming if needed), topics, 
# reps, and subsets to prep the list of model inputs, that will actually be
# appended to with the results, which are done at the top with LDA_call.
# LDA_call prints the message, and then preps the input arguments as used
# in a do.call on the function given in the control list. there are
# examples of both a generalized topicmodels LDA function wrapper (works
# with Gibbs sampling approach too!) and an "identity" 1-topic function.
# these functions allow for pre- and post-processing around the main model
# function without having to impact the general LDA_call
# select_LDA does a general measuring and then selecting, which now allows
# for arguments to be passed in!


LDA <- function(data, topics = 2, reps = 1, control = list()){
  control <- do.call("LDA_control", control)
  messageq("----Linguistic Decomposition Analyses----", control$quiet)
  LDAs <- prep_LDA_models(data = data, topics = topics, reps = reps,
                          control = control)
  nLDA <- length(LDAs)
  for (i in 1:nLDA){
    LDAs[[i]] <- LDA_call(LDA = LDAs[[i]], control = control)
  }
  selected_LDAs <- select_LDA(LDAs = LDAs, control = control)
  package_LDA(selected_LDAs = selected_LDAs, LDAs = LDAs, control = control)
}


LDA_control <- function(LDA_function = topicmodels_LDA, 
                        LDA_args = list(method = "VEM", seeded = TRUE),
                        measurer_function = AIC,
                        measurer_args = list(),
                        selector_function = which.min,
                        selector_args = list(), 
                        nsubsets = 1,
                        subset_rule = NULL,
                        soften = TRUE, 
                        quiet = FALSE){
  list(LDA_function = LDA_function,  LDA_args = LDA_args, 
       measurer_function = measurer_function, measurer_args = measurer_args, 
       selector_function = selector_function, selector_args = selector_args,
       nsubsets = nsubsets, subset_rule = subset_rule,
       soften = soften, quiet = quiet)
}

prep_LDA_models <- function(data, topics = 2, reps = 1, control = list()){
  data <- conform_data(data = data, control = control)
  subsets <- names(data)
  if(length(reps) < length(topics)){
    reps <- rep(reps, length(topics))
  }
  LDA_topics <- rep(topics, reps)
  LDA_reps <- sequence(reps)
  LDA_subsets <- rep(subsets, each = length(LDA_reps))  
  LDA_reps <- rep(LDA_reps, length(subsets))
  LDA_topics <- rep(LDA_topics, length(subsets))
  nLDA <- length(LDA_topics)
  LDAs <- vector("list", length = nLDA)
  for(i in 1:nLDA){
    LDAs[[i]] <- list(data = data[[LDA_subsets[[i]]]], 
                      data_subset = LDA_subsets[[i]],
                      topics = LDA_topics[[i]], rep = LDA_reps[[i]])
  }
  names(LDAs) <- paste0("model_", 1:nLDA)
  LDAs
}


LDA_call <- function(LDA = NULL, control = list()){
  control <- do.call("LDA_control", control)  
  LDA_msg(LDA = LDA, quiet = control$quiet)
  fun <- control$LDA_function
  args <- update_list(control$LDA_args, LDA = LDA)
  if(control$soften){
    tryCatch(do.call(what = fun, args = args), 
             warning = function(x){eval(x$call)}, 
             error = function(x = list()){list(error = x$message)})
  } else{
    do.call(what = fun, args = args)
  }
}




LDA_msg <- function(LDA, quiet = FALSE){
  subset_msg <- paste0("  data subset ", LDA$data_subset)
  topic_msg <- paste0(", ", LDA$topics, " topics")
  rep_msg <- paste0(", replicate ", LDA$rep)
  messageq(paste0(subset_msg, topic_msg, rep_msg), quiet)
}


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
    out <- list(params = list(alpha = mod@alpha, beta = mod@beta),
                document_topic_matrix = mod@gamma, 
                test_document_topic_matrix = NULL, #not yet available
                log_likelihood = mod_ll, data = data,
                topics = topics, rep = rep, data_subset = data_subset)
    class(out) <- c("LDA", "list")
    out
  } 
}

identity_LDA <- function(LDA){
  data <- LDA$data
  rep <- LDA$rep
  data_subset <- LDA$data_subset
  document_topic_table <- data$train$document_term_table 
  document_topic_table <- document_topic_table / rowSums(document_topic_table)
  colnames(document_topic_table) <- NULL
  out <- list(params = list(), document_topic_table = document_topic_table, 
              log_likelihood = NULL, data = data,
              topics = 1, rep = rep, data_subset = data_subset)
  class(out) <- c("LDA", "list")
  out
}


AIC.LDA <- function(object, ..., k = 2){
  lls <- logLik(object)
  -2 * as.numeric(lls) + k * attr(lls, "df")
}

logLik.LDA <- function(object, ...){
  object$log_likelihood
}



select_LDA <- function(LDAs = list(), control = list()){
  nLDAs <- length(LDAs)
  maxtopics <- 0
  for(i in 1:nLDAs){
    maxtopics <- max(c(maxtopics, LDAs[[i]]$topics))
  }
  if(maxtopics == 1){
    return(LDAs)
  }
  vals <- measure_LDA(LDAs = LDAs, control = control)
  fun <- control$selector_function
  args <- update_list(control$selector_args, x = vals)
  selection <- do.call(what = fun, args = args)
  LDAs[selection]  
}

measure_LDA <- function(LDAs = list(), control = list()){
  fun <- control$measurer_function
  args <- control$measurer_args
  nLDAs <- length(LDAs)
  vals <- rep(NA, nLDAs)
  for(i in 1:nLDAs){
    args <- update_list(args, object = LDAs[[i]])
    vals_i <- do.call(what = fun, args = args)
    if(length(vals_i) != 0){
      vals[i] <- vals_i
    }
  }
  vals
}

package_LDA <- function(selected_LDAs, LDAs, control = list()){
  out <- list(selected_LDAs = selected_LDAs, LDAs = LDAs, control = control)
  class(out) <- c("LDA_set", "list")
  out
}
