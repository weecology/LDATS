LDA_TS

  LDA_TS_control
    LDA_control
    TS_control

  LDA
    LDA_control
    prep_LDA_models
    LDA_call 
      - uses control$model
    select_LDA
      measure_LDA
        - uses control$measurer
      - uses control$selector
    package_LDA

  TS
    TS_control
    prep_TS_models
    sequential_TS
      sequential_TS_control
      est_changepoints
        - uses control$method
      est_regressors
        - uses TS$reponse
      package_sequential_TS
    select_TS
      measure_TS
        - uses control$measurer
      - uses control$selector
    package_TS

  package_LDA_TS



to consider
TS versus control object passing around...like once the controls are in the 
    model object why still pass them around>??
make sequential_TS the model like in LDA_call
condense x_call
bring select_x within package_x
rename prep_x_models as prep_x