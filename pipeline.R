LDA_TS

  LDA_TS_control
    LDA_control
    TS_control

  LDA
    prepare_LDA
      LDA_control
      conform_data
    run_LDA
      LDA_call (replicated for each model)
       LDA_msg
       LDA$control$model with LDA$control$model_args
    package_LDA
      measure_LDA (replicated for each model)
       LDA$control$measurer with LDA$control$measurer_args
      select_LDA
       LDA$control$selector with LDA$control$selector_args

  TS
    prepare_TS
      TS_control
    run_TS
      TS_call (replicated for each model)
       TS_msg
       TS$control$model with TS$control$model_args

[one more level deeper here!]

    package_TS
      measure_TS (replicated for each model)
       TS$control$measurer with TS$control$measurer_args
      select_TS
       TS$control$selector with TS$control$selector_args




