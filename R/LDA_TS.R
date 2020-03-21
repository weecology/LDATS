#' @title Run a set of Linguistic Decomposition Analysis models coupled to
#'   Bayesian Time Series models
#' 
#' @description The main interface function for analyzing compositional 
#'   time series using the LDATS application of Linguistic Decomposition 
#'   Analysis and Time Series modeling generally following Christensen 
#'   \emph{et al.} (2018). 
#'
#' @details For a (potentially subset) dataset consisting of counts of words 
#'   across multiple documents in a corpus, 
#'   \enumerate{
#'     \item Conduct multiple Linguistic Decomposition Analysis (LDA) models 
#'       (e.g., Latent Dirichlet Allocation using the Variational Expectation
#'       Maximization (VEM) algorithm; Blei \emph{et al.} 2003), 
#'     \item Select from the LDA model results to pick those used in the Time
#'       Series (TS) models,
#'     \item Conduct multiple compositional Bayesian TS models 
#'       (e.g., changepoint softmax regression; Ripley 1996, Venables 
#'       and Ripley 2002, Western and Kleykamp 2004, Bishop 2006, Ruggieri 
#'       2013) via a generalized linear modeling approach (McCullagh and 
#'       Nelder 1989) and using parallel tempering Markov Chain Monte Carlo
#'      (ptMCMC) methods (Earl and Deem 2005),
#'     \item Select from the TS model results to pick those used to summarize
#'       the whole model, and
#'     \item Package the results.
#'   }
#'
#' @param formulas Vector of \code{\link[stats]{formula}}(s) defining the 
#'   regression between the change points. Any predictor variable included 
#'   must also be a column in \code{data} and any (compositional) response 
#'   variable must be a set of columns in \code{data}. \cr
#'   Each element (formula) in the vector is evaluated for each number of 
#'   change points and each LDA model. \cr
#'   (See \code{\link{TS}}.)
#'
#' @param nchangepoints \code{integer}-conformable vector corresponding to the 
#'   number of change points to include in the models. 0 is valid (corresponds
#'   to no change points, so a singular time series model) and the current 
#'   implementation can reasonably include up to 6 change points. The 
#'   number of change points is used to dictate the segmentation of the 
#'   time series into chunks fit with separate models dictated by 
#'   \code{formula}. \cr
#'   Each element in the vector is the number of change points 
#'   used to segment the data for each formula (entry in \code{formulas}) 
#'   component of the TS model, for each selected LDA model. \cr
#'   (See \code{\link{TS}}.)
#'
#' @param timename \code{character} element indicating the time variable
#'   used in the time series. Defaults to \code{"time"}. The variable must be
#'   integer-conformable or a \code{Date}. If the variable named
#'   is a \code{Date}, the input is converted to an integer, resulting in the
#'   timestep being 1 day, which is often not desired behavior. \cr
#'   (See \code{\link{TS}}.)
#'
#' @param weights Optional class \code{numeric} vector of weights for each 
#'   document. Defaults to \code{NULL}, translating to an equal weight for
#'   each document. When using \code{\link{sequential_TS}} in a standard LDATS 
#'   analysis, it is advisable to weight the documents by their total size,
#'   as the result of, e.g., \code{\link[topicmodels]{LDA}} is a matrix of 
#'   proportions, which does not account for size differences among documents.
#'   For most models, a scaling of the weights (so that the average is 1) is
#'   most appropriate, and this is accomplished using \code{document_weights}.
#'   \cr
#'   (See \code{\link{TS}}.)
#'
#' @param control \code{list} of parameters to control the fitting of the
#'   LDATS model. Values not input assume defaults set by 
#'   \code{\link{LDA_TS_control}}.
#'
#' @param data Any of the data structures allowable for LDATS analyses:
#'   \code{matrix} or \code{data.frame} document term table, 
#'   \code{list} of document term and covariate tables, a \code{list} of 
#'   training and test sets of the two tables, or a \code{list} of multiple 
#'   replicate splits of training and test sets of the two tables. \cr
#'   See \code{\link{conform_data}}, which is used to ensure data structure
#'   validity for the desired model.
#'  
#' @param topics \code{integer}-conformable \code{vector} of the number of 
#'   topics to evaluate for each model. \cr
#'   (See \code{\link{LDA}}.)
#'
#' @param reps \code{integer}-conformable number of replicate starts to use 
#'   for each value of \code{topics}. \cr
#'   (See \code{\link{LDA}}.)
#'
#' @return \code{LDA_TS} \code{list} with all fitted LDA and TS models and 
#'   selected models specifically as elements named
#'   \describe{
#'     \item{\code{LDA models}}{\code{list} of all and selected models as well
#'       as controls from \code{\link{LDA}}}
#'     \item{\code{TS models}}{\code{list} of all and selected models as well
#'       as controls from \code{\link{TS}}}
#'     \item{\code{control}}{\code{list} of overall model controls}
#'   }
#' 
#' @references 
#'   Blei, D. M., A. Y. Ng, and M. I. Jordan. 2003. Latent Dirichlet
#'   Allocation. \emph{Journal of Machine Learning Research} 
#'   \strong{3}:993-1022.
#'   \href{http://jmlr.csail.mit.edu/papers/v3/blei03a.html}{link}.
#'
#'   Bishop, C. M. 2006. \emph{Pattern Recognition and Machine Learning}. 
#'    Springer, New York, NY, USA.
#'
#'   Christensen, E., D. J. Harris, and S. K. M. Ernest. 2018.
#'   Long-term community change through multiple rapid transitions in a 
#'   desert rodent community. \emph{Ecology} \strong{99}:1523-1529. 
#'   \href{https://doi.org/10.1002/ecy.2373}{link}.
#'
#'   Earl, D. J. and M. W. Deem. 2005. Parallel tempering: theory, 
#'   applications, and new perspectives. \emph{Physical Chemistry Chemical 
#'   Physics} \strong{7}: 3910-3916.
#'   \href{https://doi.org/10.1039/B509983H}{link}.
#'
#'   Grun B. and K. Hornik. 2011. topicmodels: An R Package for Fitting Topic
#'   Models. \emph{Journal of Statistical Software} \strong{40}:13.
#'   \href{https://www.jstatsoft.org/article/view/v040i13}{link}.
#'
#'   McCullagh, P. and J. A. Nelder. 1989. \emph{Generalized Linear Models}.
#'   2nd Edition. Chapman and Hall, New York, NY, USA.
#'
#'   Ripley, B. D. 1996. \emph{Pattern Recognition and Neural Networks}. 
#'   Cambridge University Press, Cambridge, UK.
#'
#'   Ruggieri, E. 2013. A Bayesian approach to detecting change points in 
#'   climactic records. \emph{International Journal of Climatology}
#'   \strong{33}:520-528.
#'   \href{https://doi.org/10.1002/joc.3447}{link}.
#'
#'   Venables, W. N. and B. D. Ripley. 2002. \emph{Modern and Applied
#'   Statistics with S}. Fourth Edition. Springer, New York, NY, USA.
#'
#'   Western, B. and M. Kleykamp. 2004. A Bayesian change point model for 
#'   historical time series analysis. \emph{Political Analysis}
#'   \strong{12}:354-374.
#'   \href{https://doi.org/10.1093/pan/mph023}{link}.
#'
#' @export
#'
LDA_TS <- function(data, topics = 2, reps = 1, formulas = ~ 1, 
                   nchangepoints = 0, timename = "time", weights = TRUE, 
                   control = list()){
  control <- do.call("LDA_TS_control", control)
  data <- conform_data(data = data, control = control)
  LDAs <- LDA(data = data, topics = topics, reps = reps, 
              control = control$LDA_control)
  TSs <- TS(LDAs = LDAs, data = data, formulas = formulas, 
            nchangepoints = nchangepoints, timename = timename, 
            weights = weights, control = control$TS_control) 
  package_LDA_TS(LDAs = LDAs, TSs = TSs, control = control)
}


#' @title Package LDA and TS model outputs
#'
#' @description Combine the results from each model component.
#'
#' @param LDAs \code{LDA_set} \code{list} of selected and all LDAs from 
#'   \code{\link{LDA}}.
#'
#' @param TSs \code{TS_set} \code{list} of selected and all TSs from 
#'   \code{\link{TS}}.
#'
#' @param control A \code{list} of parameters to control the fitting of the
#'   LDATS model. Values not input assume defaults set by 
#'   \code{\link{LDA_TS_control}}.
#'
#' @return \code{LDA_TS} \code{list} with all fitted LDA and TS models and 
#'     selected models specifically as elements named
#'     \code{"LDA models"} (from \code{\link{LDA_set}}),
#'     \code{"Selected LDA models"} (from \code{\link{select_LDA}}), 
#'     \code{"TS models"} (from \code{\link{TS_on_LDA}}), and
#'     \code{"Selected TS models"} (from \code{\link{select_TS}}) elements.
#'
#' @export
#'
package_LDA_TS <- function(LDAs, TSs, control){
  out <- list("LDA models" = LDAs, "TS models" = TSs, control = control)
  class(out) <- c("LDA_TS", "list")
}



#' @title Create the controls list for the complete Linguistic Decomposition
#'   Analysis and Time Series model
#'
#' @description Defines and creates a \code{list} used to control the 
#'   decomposition and time series model fitting occurring within 
#'   \code{\link{LDA_TS}}. 
#'
#' @param LDA_model Main LDA \code{function}.
#'
#' @param LDA_model_args \code{list} of (named) arguments to be used in 
#'   \code{LDA_model} via \code{\link{LDA_call}}.
#'
#' @param TS_model Main TS \code{function}.
#'
#' @param TS_model_args \code{list} of (named) arguments to be used in 
#'   \code{TS_model}.
#'
#' @param TS_response \code{function} used to model the compositional 
#'   response.
#'
#' @param response_args \code{list} of (named) arguments to be used in 
#'   \code{TS_response} via \code{\link{do.call}}. 
#'   \cr \cr
#'   Could be managed via a \code{<reponse>_TS_control} function like
#'   \code{\link{multinom_TS_control}}.
#'
#' @param TS_method \code{function} used to drive the sampler of the TS
#'   models; \code{TS_method} defines and operates the computational 
#'   procedure. \cr \cr
#'   Current pre-built options include \code{\link{ldats_classic}}.
#'
#' @param TS_method_args \code{list} of (named) arguments to be used in 
#'   \code{TS_method} via \code{\link{do.call}}. 
#'   \cr \cr
#'   Could be managed via a \code{<method>_control} function like
#'   \code{\link{ldats_classic_control}}.
#'
#' @param summary_prob Probability used for summarizing the posterior 
#'   distributions (via the highest posterior density interval, see
#'   \code{\link[coda]{HPDinterval}}).
#'
#' @param quiet \code{logical} indicator of whether the model should run 
#'   quietly (if \code{FALSE}, a progress bar and notifications are printed).
#'
#' @param soften \code{logical} indicator of whether the model should error 
#'   softly or if errors should trigger a full-stop to the pipeline.
#'
#' @param TS_measurer \code{function} used in evaluation of the TS
#'   models; \code{measurer} creates a value for each model.
#'
#' @param TS_measurer_args \code{list} of (named) arguments to be used in 
#'   \code{TS_measurer} via \code{\link{do.call}}. 
#'
#' @param TS_selector \code{function} usde in evaluation of the TS
#'   models; \code{TS_selector} operates on the values to choose the models. 
#'
#' @param TS_selector_args \code{list} of (named) arguments to be used in 
#'   \code{TS_selector} via \code{\link{do.call}}. 
#'
#' @param LDA_measurer \code{function} used in evaluation of the LDA
#'   models; \code{LDA_measurer} creates a value for each model.
#'
#' @param LDA_measurer_args \code{list} of (named) arguments to be used in 
#'   \code{LDA_measurer} via \code{\link{do.call}}. 
#'
#' @param LDA_selector \code{function} usde in evaluation of the LDA
#'   models; \code{LDA_selector} operates on the values to choose the models. 
#'
#' @param LDA_selector_args \code{list} of (named) arguments to be used in 
#'   \code{LDA_selector} via \code{\link{do.call}}. 
#'
#' @param ... Not passed along to the output, rather included to allow for
#'   automated removal of unneeded controls.
#'
#' @param nsubsets Number of data subsets.
#'
#' @param subset_rule \code{function} used to subset the data.
#'
#' @return \code{list} of \code{list}s and single elements, with named 
#'   elements corresponding to the arguments.
#'
#' @examples
#'   LDA_TS_control()
#'
#' @export
#'
LDA_TS_control <- function(LDA_model = topicmodels_LDA, 
                           LDA_model_args = 
                             list(method = "VEM", seeded = TRUE),
                           LDA_measurer = AIC,
                           LDA_measurer_args = list(),
                           LDA_selector = which.min,
                           LDA_selector_args = list(), 
                           TS_model = sequential_TS,
                           TS_model_args = sequential_TS_control(),
                           TS_response = multinom_TS,
                           TS_response_args = multinom_TS_control(),
                           TS_method = ldats_classic,
                           TS_method_args = ldats_classic_control(),
                           nsubsets = 1,
                           subset_rule = NULL,
                           summary_prob = 0.95,
                           soften = TRUE, 
                           quiet = FALSE,
                           TS_measurer = AIC,
                           TS_measurer_args = list(),
                           TS_selector = which.min,
                           TS_selector_args = list(), ...){

  LDA_control <- LDA_control(model = LDA_model, model_args = LDA_model_args, 
                             measurer = LDA_measurer, 
                             measurer_args = LDA_measurer_args, 
                             selector = LDA_selector, 
                             selector_args = LDA_selector_args,
                             nsubsets = nsubsets, subset_rule = subset_rule,
                             soften = soften, quiet = quiet)

  TS_control <- TS_control(response = TS_response, 
                           response_args = TS_response_args, 
                           method = TS_method, 
                           method_args = TS_method_args, 
                           measurer = TS_measurer,
                           measurer_args = TS_measurer_args, 
                           selector = TS_selector, 
                           selector_args = TS_selector_args,
                           summary_prob = summary_prob, 
                           soften = soften, quiet = quiet)

  list(LDA_control = LDA_control, TS_control = TS_control, 
       nsubsets = nsubsets, subset_rule = subset_rule,
       soften = soften, quiet = quiet)
}


