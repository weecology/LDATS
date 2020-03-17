time_order_data ... hmmmmm


standardized output from LDA functions
s3 object with
params (however they're provided, collated into a list)
document_topic_table
log_likelihood as a logLik object
data used to fit the model

def still in the thick of this
just need to take a solid break 
currently working in LDA to get it to work with function args


LDA_set is now LDA and TS_on_LDA is now TS

trying to really streamline the code at each level of the pipeline and such

now LDA could be any of the linquistic decomposition analyses or whatever
including any of a number of "LDA" functions or models, so i'm revoking
the importing of LDA from topicmodels but keeping topicmodels imported
to allow calling of topicmodels::LDA from inside LDA (the default!)


generalizing multinom_TS to be compositional_TS

in the multi_LDA_TS
  in test_LDA_TS
  requires a predict.LDA_TS

LDAT_TS
  developing predict method

check seeds is now check_nreps
seed is rep

introduction of soften logical variable
designed to help soften errors in pipelines through wrapping in tryCatch