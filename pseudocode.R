

# pipeline

meta_LDA_TS
  for each iteration
   - split_LDA_TS_data (divides data into training and test sets)
   - LDA_TS (runs on the training set)
   - test_LDA_TS (measures the model on the test data)
  measure_LDA_TS  (combines all the individual iterations' measures)
  select_LDA_TS (only used if multiple models are given)
  package_meta_LDA_TS

run meta_LDA_TS for each model
iterates over the subsets of the data for a given model
(or you could give it multiple models and it will run using all of them...
  for this you'll want to be able to pass multiple models between the stages)

