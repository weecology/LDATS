# still needs to be evolved and adapted with usages!


# if the data are a data.frame or matrix, or list of them, or list of lists
# of them, you're good. output is a list of lists of data.frames
# this allows there to be iterations of the same data set, etc.
# or actually one more level of the list because of the data splitting lololz
# data <- rodents
# conform_data(data)
# 
# conform_data(data)
conform_data <- function(data, control = list()){

  depth <- list_depth(data)
  if(depth == 0){  
    if(inherits(data, "data.frame") | inherits(data, "matrix")){
      msg <- "covariate table not provided, assuming equi-spaced data"
      messageq(msg, control$quiet)
      nobs <- nrow(data)
      covariate <- data.frame(time = 1:nobs)
      data <- list(document_term_table = data, 
                   document_covariate_table = covariate)
      depth <- list_depth(data)
    } else{
      stop("improper data format")
    }
  }
  if(depth == 1){
    which_term <- grep("term", names(data), ignore.case = TRUE)
    which_covariate <- grep("covariate", names(data), ignore.case = TRUE)
    if(length(which_term) != 1){
      stop("one, and only one, element in `data` can include `term`")
    }
    if (length(which_covariate) == 0){
      msg <- "covariate table not provided, assuming equi-spaced data"
      messageq(msg, control$quiet)
      nobs <- nrow(data[[which_term]])
      covariate <- data.frame(time = 1:nobs)
      data$document_covariate_table <- covariate    
    } else if(length(which_covariate) > 1){
      stop("at most one element in `data` can include `covariate`")
    }
    names(data)[which_term] <- "document_term_table"
    names(data)[which_covariate] <- "document_covariate_table"
    data <- list(data)
    depth <- list_depth(data)
  } 
  if(depth == 2){
    nsubsets_in <- length(data)
    nsubsets_out <- control$nsubsets

    if(nsubsets_in != 1 && nsubsets_in != nsubsets_out){
      stop("mimatched request for data subsets")
    }

    if(nsubsets_out > 0){  
      data_1 <- data[[1]]
      for(i in 1:nsubsets_out){

        rule <- control$rule
        if(is.null(rule)){
          rule <- null_rule
        }
        if(nsubsets_in == 1){
          data_i <- data_1
        } else if (nsubsets_in > 1){
          data_i <- data[[i]]
        }


        if(!all(c("test", "train") %in% names(data_i))){
          which_term <- grep("term", names(data_i), ignore.case = TRUE)
          which_covariate <- grep("covariate", names(data_i), 
                                  ignore.case = TRUE)
          if(length(which_term) != 1){
            stop("one, and only one, element in `data` can include `term`")
          }
          if (length(which_covariate) == 0){
            msg <- "covariate table not provided, assuming equi-spaced data"
            messageq(msg, control$quiet)
            nobs <- nrow(data_i[[which_term]])
            covariate <- data.frame(time = 1:nobs)
            data_i$document_covariate_table <- covariate    
          } else if(length(which_covariate) > 1){
            stop("at most one element in `data` can include `covariate`")
          }
          names(data_i)[which_term] <- "document_term_table"
          names(data_i)[which_covariate] <- "document_covariate_table"
        }
        dtt <- data_i$document_term_table
        dct <- data_i$document_covariate_table
        args <- list(data = dtt, iteration = i)
        test_train <- do.call(what = rule, args = args)
        in_train <- test_train == "train"
        in_test <- test_train == "test"
        train <- list(document_term_table = dtt[in_train, ],
                      document_covariate_table = dct[in_train, ])
        test <- list(document_term_table = dtt[in_test, ],
                     document_covariate_table = dct[in_test, ])
        data[[i]] <- list(test = test, train = train)
      }
      names(data) <- paste0("subset_", 1:nsubsets_out)
    }
    depth <- list_depth(data)
  }
  if(depth == 3){
    
  }
  data

}
