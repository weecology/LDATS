#' @title Conform data for LDATS modeling
#'
#' @description Given any of a variety of possible data input types 
#'  (\code{data.frame}/\code{matrix}, \code{list}, \code{list} of 
#'  \code{list}s, or \code{list} of \code{list} of \code{list}s) and 
#'  controls, this produces a properly formatted set of data (sets) for 
#'  LDATS modeling.
#'
#' @details This function makes use of the \code{\link{list_depth}}
#'  utility that recursively works through an object to tell you
#'  how nested a lists is. Working up from the most elemental version
#'  possible, if it's not a \code{list}, but the data are a term table, the
#'  covariate table is added with assumed equispersed data like before and
#'  the data are now a \code{list}. Then, if it is a \code{list} but only a 
#'  of depth 1 (a \code{list} of two tables), we need to wrap it in a 
#'  \code{list} to make it depth-2, functionally a 1-subset data set. Then, 
#'  for a \code{list} of depth two, we need to potentially expand to a 
#'  multiple-subset data set, to allow for cross validtion methods, for 
#'  example. So, the \code{list} of depth 2 is replicated out to create a 
#'  longer \code{list} that is still depth 2 but is now of length 
#'  \code{control$nsubsets}. Then, the subsetting of the data occurs 
#'  according to the \code{control$subset_rule}, and each depth-2 \code{list}
#'  is actually split to a final level of train and test subsets of the data,
#'  making the \code{list} depth 3. \cr \cr
#'  The training and testing data are saved as trimmed versions of the 
#'  two tables. 
#'
#' @param data A document term table, \code{list} of document term and 
#'   covariate tables, a list of training and test sets of the two tables,
#'   or a list of multiple replicate splits of training and test sets of
#'   the two tables. 
#'
#' @param control \code{list} of control options for the data conforming.
#'
#' @return \code{list} of properly formatted LDATS data.
#'
#' @export
#'
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

#' @title Subset data sets 
#'
#' @description For use within, e.g., cross validation methods, these 
#'  functions subdivide the data into testing and training subsets. \cr \cr
#'  \code{null_rule} places all data in the training set. \cr \cr 
#'  \code{random_loo} conducts randomized leave-one-out with no buffer. 
#'   \cr \cr 
#'  \code{systematic_loo} conducts systematic leave-one-out with no buffer. 
#'   Assumes 1:1 between iteration and datum location to drop. \cr \cr 
#'  \code{leave_p_out} is a fully flexible leave p out function allowing for 
#'   asymmetric buffers and randomization. If \code{random = TRUE}, the test 
#'   data are selected randomly, otherwise locations are used.
#'
#' @param data \code{data.frame} or \code{matrix} of data to be split.
#'
#' @param iteration \code{integer}-conformable value indicating which 
#'  iteration through the process the current implementation is.
#'
#' @param p \code{integer}-conformable value of how many samples to leave out.
#'
#' @param pre,post \code{integer}-conformable values of how many samples
#'  to include in the buffer around the focal left out data. Can be 
#'  asymmetric.
#' 
#' @param random \code{logical} indicator of if the left out data should be
#'  randomly selected.
#'
#' @param locations \code{integer}-conformable values referencing which
#'  data to hold out.
#'
#' @return \code{character} \code{vector} of \code{"train"} and \code{"test"}
#'  values.
#'
#' @name data_subsetting
#'


#' @rdname data_subsetting
#'
#' @export
#'
null_rule <- function(data, iteration = 1){
  n <- NROW(data)
  rep("train", n)
}


#' @rdname data_subsetting
#'
#' @export
#'
systematic_loo <- function(data, iteration = 1){
  leave_p_out(data = data, random = FALSE, locations = iteration)
}

#' @rdname data_subsetting
#'
#' @export
#'
random_loo <- function(data, iteration = 1){
  leave_p_out(data = data)
}

#' @rdname data_subsetting
#'
#' @export
#'
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

