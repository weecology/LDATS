#' Run a batch of LDA models in parallel
#' 
#' Runs each of the number of input topics for the number of seeds
#' 
#' @param data Table of integer data (species counts by period)
#' @param ntopics set of topics to evaluate
#' @param nseeds number of seeds to use
#' @param method what LDA method to use (currently only support for "VEM")
#' @param sort Logical of whether or not to sort the output
#' @param sortby Which metric to sortby ("aic" or "aicc")
#' @param parallel Logical of whether or not to run the models in parallel
#' @param ncores Integer number of cores to use
#' 
#' @return List of [1] list of models and [2] table model summaries
#' 
#' @examples aic_values = batch_LDA(data, 2, 1, "VEM", TRUE, "aicc", FALSE)
#' @export 


batch_LDA <- function(data, ntopics, nseeds, method, sort, sortby,
                            parallel, ncores = NULL){

  # verify method, currently only supporting VEM

    if(method != "VEM"){
      return(paste("Method ", method, " not currently supported"))
    }

  # set up the seeds and numbers of topics to cycle through

    seed_in <- rep(seq(2, nseeds * 2, 2), length(ntopics))
    k_in <- rep(ntopics, each = length(seq(2, nseeds * 2, 2)))

 
  # if parallel, run in parallel

    if(parallel == TRUE){

      if(length(ncores) == 0){
        ncores <- floor(0.75 * detectCores(logical = FALSE))
      }

      # register the parallel cluster

        cl <- parallel::makeCluster(ncores)
        doParallel::registerDoParallel(cl)

      # run the models

        models <- foreach::foreach(i = 1:length(seed_in), 
                                   .packages = "topicmodels") %dopar% {
          topicmodels::LDA(data, k = k_in[i], 
                           control = list(seed = seed_in[i]),
                           method = 'VEM')
        }

        names(models) <- paste("seed_", seed_in, "_ntopics_", k_in, sep = "")

      # summarize the models

        modelsummaries <- foreach::foreach(i = 1:length(seed_in), 
                                           .packages = "topicmodels",
                                           .combine = rbind) %dopar% {
          x <- models[[i]]
          llik <- sum(x@loglikelihood)
          postx <- modeltools::posterior(x) 
          nparam <- (nrow(postx$topics) - 1) * (ncol(postx$topics)) + 
                     nrow(postx$terms) * (ncol(postx$terms) - 1) 
          nobs <- nrow(data) * ncol(data)
          aic <- 2 * nparam - 2 * llik
          aicc <- aic + 
                   (((2 * (nparam ^ 2)) + (2 * nparam)) / (nobs - nparam - 1))

          c(seed = seed_in[i], k = k_in[i], llik = llik,
            nobs = nobs, nparam = nparam, aic = aic, aicc = aicc)

        }
  
      # if wanted, sort the output

        if(sort == TRUE){

          # verify the metric for sorting and sort

            if(sortby %in% c("aic", "aicc")){
              modelsT <- models
              modelsummariesT <- modelsummaries
              om <- order(modelsummariesT[,sortby])
              modelsummaries <- modelsummariesT[om, ]
              rownames(modelsummaries) <- om

              models <- foreach::foreach(i = 1:length(seed_in), 
                                         .packages = "topicmodels") %dopar% {
                modelsT[[om[i]]]
              }
              names(models) <- names(modelsT)[om]
            }else{

              print(paste("Sortby input ", sortby, 
                          " not recognized. Output not sorted", sep = ""))
            }

        }  

      parallel::stopCluster(cl)

    }else{

      # prep the output

        models <- vector("list", length(seed_in))
        names(models) <- paste("seed_", seed_in, "_ntopics_", k_in, sep = "")
        modelsummaries <- data.frame(seed = seed_in,
                                     k = k_in,
                                     llik = rep(NA, length(seed_in)),
                                     nobs = rep(NA, length(seed_in)),
                                     nparam = rep(NA, length(seed_in)),  
                                     aic = rep(NA, length(seed_in)),
                                     aicc = rep(NA, length(seed_in)))


      # iterate through, run and summarize

        for(i in 1:length(seed_in)){

          models[[i]] <- topicmodels::LDA(data, k = k_in[i], 
                                          control = list(seed = seed_in[i]),
                                          method = 'VEM')

          x <- models[[i]]
          llik <- sum(x@loglikelihood)
          postx <- modeltools::posterior(x) 
          nparam <- (nrow(postx$topics) - 1) * (ncol(postx$topics)) + 
                     nrow(postx$terms) * (ncol(postx$terms) - 1) 
          nobs <- nrow(data) * ncol(data)
          aic <- 2 * nparam - 2 * llik
          aicc <- aic + 
                   (((2 * (nparam ^ 2)) + (2 * nparam)) / (nobs - nparam - 1))

          modelsummaries[i, ] <- c(seed = seed_in[i], k = k_in[i], 
                                   llik = llik, nobs = nobs, nparam = nparam,
                                   aic = aic, aicc = aicc)

        }

      # if wanted, sort the output

        if(sort == TRUE){

          # verify the metric for sorting and sort

            if(sortby %in% c("aic", "aicc")){
              modelsT <- models
              modelsummariesT <- modelsummaries
              om <- order(modelsummariesT[,sortby])
              modelsummaries <- modelsummariesT[om, ]
              rownames(modelsummaries) <- om

              for(i in 1:length(seed_in)){
                models[[i]] <- modelsT[[om[i]]]
              }
              names(models) <- names(modelsT)[om]

            }else{

              print(paste("Sortby input ", sortby, 
                          " not recognized. Output not sorted", sep = ""))
            }

        }  

    }


  # prep output 

    output <- list(models, modelsummaries)
    names(output) <- c("Models", "ModelSummaries")

  # return output

    return(output)

}
