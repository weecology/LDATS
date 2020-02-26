# data can come into LDA_TS LDA TS in a variety of forms, and depending on 
# usages, might take a variety of different forms
# the purpose of this function is to generalize and extract the code used
# to shuddle between data formats from functions / replace with a single line
# it's still a work in progress and needs more extensive usage exploration,
# as it's going to be a workhorse function. 
# this function makes use of a utility i brought over from portalcasting
# called list_depth that recursively works through an object to tell you
# how nested lists are. its extremely useful when you could have a list or
# a list of multiple lists and need to easily distinguish
# the idea is as follows: working up from the most elemental version
# possible, if it's not a list, but the data are a term table, the covariate
# table is added with assumed equispersed data like before and the data are
# now a list
# then, if it is a list but only a list of depth 1 (a list of two tables)
# we actually need to wrap it in a list to make it depth 2...think of this
# as a 1-subset data set. then, for a list of depth two, we need to 
# potentially expand to a multiple-subset data set, to allow for cross valid
# methods, for example. so the list of depth 2 is replicated out to 
# create a longer list that is still depth 2 but is now of length 
# control$nsubsets. and then the subsetting of the data occurs according to 
# the control$subset_rule, and each depth-2 list is actually split to
# a final level of train and test data, making the list depth 3
# the training and testing data are saved as trimmed versions of the 
# two tables. currently its not saving the test/train split explicitly,
# just implicitly via the data encoding that exists. we should probably
# shore this up a bit more for sure.
# also this function is big and modularized a good degree already...it could
# get chunked into subfunctions
# there are functions for basic leave p out cross validation, including
# both systematic and random approaches
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
      names(data) <- 1:nsubsets_out
    }
    depth <- list_depth(data)
  }
  if(depth == 3){

    nsubsets_in <- length(data)
    nsubsets_out <- control$nsubsets

    if(nsubsets_in != 1 && nsubsets_in != nsubsets_out){
      stop("mimatched request for data subsets")
    }

    for(i in 1:nsubsets_out){

      data_i <- data[[i]]
      which_train <- grep("train", names(data_i), ignore.case = TRUE)
      which_test <- grep("test", names(data_i), ignore.case = TRUE)
      if(length(which_train) != 1){
        stop("only one, element in a `data` subset can include `train`")
      }
      if(length(which_test) != 1){
        stop("only one, element in a `data` subset can include `test`")
      }
      for(j in 1:2){
        data_ij <- data[[i]][[j]]
        which_term <- grep("term", names(data_ij), ignore.case = TRUE)
        which_covariate <- grep("covariate", names(data_ij), 
                                ignore.case = TRUE)
        if(length(which_term) != 1){
          stop("one, and only one, element in `data` can include `term`")
        }
        if (length(which_covariate) == 0){
          stop("covariate table not provided, can't be made from split data")
        } else if(length(which_covariate) > 1){
          stop("at most one element in `data` can include `covariate`")
        }
      }
    }
  }
  data

}


# data subsetting rules

null_rule <- function(data, iteration = 1){
  n <- NROW(data)
  rep("train", n)
}

# simple leave one outs with no buffer

# assumes 1:1 between iteration and datum location to drop

systematic_loo <- function(data, iteration = 1){
  leave_p_out(data = data, random = FALSE, locations = iteration)
}

# randomly selected 

random_loo <- function(data, iteration = 1){
  leave_p_out(data = data)
}

# fully flexible leave p out function allowing for buffers
# if random the test data are selected randomly, otherwise locations are used

leave_p_out <- function(data, p = 1, pre = 0, post = 0, 
                        random = TRUE, locations = NULL){
  n <- NROW(data)
  test_train <- rep("train", n)

  if(random){
    locations <- sample(1:n, p)
  }

  for(i in 1:p){
    hold_out <- (locations[i] - pre):(locations[i] + post)
    test_train[hold_out] <- "out"
  }
  test_train[locations] <- "test"
  test_train
}

