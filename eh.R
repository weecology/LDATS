
select_LDA <- function(LDA_models = NULL, control = list()){
  if("LDA_set" %in% attr(LDA_models, "class") == FALSE){
    stop("LDA_models must be of class LDA_set")
  }
  control <- do.call("LDA_set_control", control)
  measurer <- control$measurer
  selector <- control$selector  
  lda_measured <- vapply(LDA_models, measurer, 0) %>%
                  matrix(ncol = 1)
  lda_selected <- apply(lda_measured, 2, selector) 
  which_selected <- which(lda_measured %in% lda_selected)
  out <- LDA_models[which_selected]
  class(out)  <- c("LDA_set", "list") 
  out
}
