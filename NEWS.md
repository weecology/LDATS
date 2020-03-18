# LDATS (development version)

Version numbers follow [Semantic Versioning](https://semver.org/).


# LDATS 0.2.7(https://github.com/weecology/ldats/releases/tag/v0.2.7)
*2020-03-18* 

## Patching CRAN issues with vignette building
* Dependencies are being managed different now
* For the paper comparison vignette, all of the code is pre-run and saved in the LDATS-replications repository
* Allows removal of otherwise unused packages from this package's dependency list

# LDATS 0.2.6(https://github.com/weecology/ldats/releases/tag/v0.2.6)
*2020-03-02* 

## Patching a bug in tests for r-devel
* `straingsAsFactors` update
* only involved patching one test

# LDATS 0.2.5(https://github.com/weecology/ldats/releases/tag/v0.2.5)
*2019-12-22* 

## General editing of simulation functions
* Don't need to make a sparse matrix to pass in now
* Tweaking the simulation functions to simplify X

## Patching a bug in `sim_TS`
* Using only as.matrix() fails if there is only 1 year in a segment and there are multiple covariates. In that case, as.matrix(X[in1, ]) returns a matrix of n_covariates rows x 1 column, instead of a matrix of 1 row and n_covariates columns. This edit should fix that by forcing it into a matrix of the correct number of rows.

# [LDATS 0.2.4](https://github.com/weecology/ldats/releases/tag/v0.2.4)
*2019-07-28*

## Edits for resubmission to CRAN
* Added space between braces for auto-linking in the description text.

# [LDATS 0.2.3](https://github.com/weecology/ldats/releases/tag/v0.2.3)
*2019-07-24*

## Edits for resubmission to CRAN
* Output given by `print`/`cat` has been replaced with `message` messages.
* Added examples in documentation (and replacement of `duntrun{}` with `donttest{}`)
* Editing of description file for specs
* Reduction of test runtimes

## Changed functions
* `messageq` replaces `qprint`

# [LDATS 0.2.2](https://github.com/weecology/ldats/releases/tag/v0.2.2)
*2019-07-10*

## Edits for submission to CRAN
* Including appropriate files in the Rbuildignore

## Minor patching of vignette code
* Handling the downloads so they work robustly locally

# [LDATS 0.2.1](https://github.com/weecology/ldats/releases/tag/v0.2.1)
*2019-07-09*

## Vignette update
* Incorporates Hao's feedback and edits on the paper comparison vignette
* Updates the vignette to work with the contemporary version of the package
* Allowed removal of the large model cache files

## Zenodo json
* Inclusion of the json file for the Zenodo page

## Tidying of the model doc
* The .pdf describing the model (the manuscript work in progress) is now at the top level and named "LDATS_model.pdf", to allow the full model description to remain stable while the ms development happens elsewhere.

# [LDATS 0.2.0](https://github.com/weecology/ldats/releases/tag/v0.2.0)
*2019-07-09*

## API updates
* At the `LDA_TS` function level, the separate inputs for data tables (`document_term_table` and `document_covariate_table`) have been merged into a single input `data`, which can be just the `document_term_table` or a list including the `document_term_table` and optionally also a `document_covariate_table`. If covariates aren't provided, the function now constructs a covariate table assuming equi-spaced observations. If using a list, the function assumes that one and only one element of the list will have a name containing the letters "term", and at most one element containing the letters "covariate" (regular expressions are used for matching). ([addresses issue 119](https://github.com/weecology/LDATS/issues/119))
* `timename` has been moved from within the `TS_controls_list` to a main argument in all associated functions.
* The control lists have been made easier to interact with. Primarily, the arguments that previously required `LDA_controls_list`, `TS_controls_list`, or `LDA_TS_controls_list` inputs now take general `list` inputs (so  `LDA_TS` does not need to have a nested set of control functions). Each control list is passed through a function (`LDA_set_control`, `TS_control`, or `LDA_TS_control`) to set any non-input values to their defaults. This also allows the removal of those controls list class definitions. ([addresses issue 130](https://github.com/weecology/LDATS/issues/130))

## Fixed and updated example code to improve user experience
* Reduced the complexity of the example in the README ([addresses issue 115](https://github.com/weecology/LDATS/issues/115))
* Added `control` input in the `plot` call in the example in the README ([addresses issue 116](https://github.com/weecology/LDATS/issues/116))
* Reduced the number of seeds in the rodent vignette example ([addresses issue 117](https://github.com/weecology/LDATS/issues/117))

## Updated calculation of the number of observations in LDA
* The number of observations for a VEM-fit LDA is now calculated as the number of entries in the document-term matrix (following Hoffman et al. and Buntine, see `?logLik.LDA_VEM` for references.
* Associated, we now include an AICc function that is general and works in this specific case as defined  ([addresses issue 129](https://github.com/weecology/LDATS/issues/129))

## Fixed bug in plotting across multiple outputs
* A few plotting functions use `devAskNewPage` to help flip through multiple outputs, but were only resetting it with `devAskNewPage(FALSE)` at the end of a clean execution. The code has been updated with `on.exit(devAskNewPage(FALSE))`, which accounts for failed executions. ([addresses issue 118](https://github.com/weecology/LDATS/issues/118))

## Renamed functions
* `summarize_TS` has been renamed `package_TS` to align with the other `package_` functions that package model output.

## Simulate functions
* Basic simulation functionality has been added for help with generating data sets to analyze. ([addresses issue 114](https://github.com/weecology/LDATS/issues/114))
* `sim_LDA_data` simulates an LDA model's document-term-matrix
* `sim_TS_data` simulates an TS model's document-topic distribution matrix
* `sim_LDA_TS_data` simulates an LDA_TS model's document-term-matrix
* `softmax` and `logsumexp` are added as utility functions

## Improved pkgdown site
* Function organization ([addresses issue 122](https://github.com/weecology/LDATS/issues/122)) and navbar formatting.

## Editing of output from `TS`
* Due to a misread of earlier code, the AIC value in the output from `TS` was named "deviance". The output has been updated to return the AIC.

## Replacement of `AIC` method with `logLik` method for `TS_fit`

# [LDATS 0.1.0](https://github.com/weecology/LDATS/pull/105)
*2019-02-11*

## Code structure
* Creation of a [standard API and code pipeline](https://weecology.github.io/LDATS/articles/LDATS_codebase.html) for all components of the LDATS analysis.
* Substantial refactor of the underlying code from hardcoded to generalized functions.
* Development of checking functions used to run the basic structural checks on the function inputs.
* Inclusion of control options lists for the LDA stage, TS stage, and overall to reduce the length of input lists.

## Full inclusion of functions
* All functions used in the code base are now exported, documented, and tested.

## LDA model AIC calculation
* `AIC.LDA_VEM()` now uses the number of parameters as reported from `logLik` to calculate AIC.
* Previous by-hand calculations of AIC included variational parameters that are integrated out of the model in the total parameter count.

## Regressor estimates
* Time series models allow for flexible covariate set for regression via formula inputs to the top-level functions.
* The time series model code now also includes estimation of the parameters defining the between-change point regressions (*i.e.*, the regressors).
* Regressor estimates come as marginal posterior distributions, and are calculated by unconditioning the estimates generated under known change points. 

## Document weighting
* `document_weights()` function is provided to allow for appropriate weighting of documents by their sizes (number of words) so that an average-length document is 1.
* Document weighting is done automatically by default, which is easily undone by using `weights = NULL`. 

## ptMCMC functionality
* The ptMCMC code has been refactored into functions, many of which are generalized to use in other contexts.
* The temperature schema is fully controllable via arguments to the TS controls list
* Burn-in removal and thinning of final chains is controllable via the TS controls list

## Optional memoisation
* Memoisation of `multinom_TS()` and `multinom_TS_chunk()` now is optional via `memoise_fun()` and is controlled through the TS controls list.

## Plotting functions
* `LDA_set()`, `LDA_TS()`, and `TS()` now all have default plotting options on their outputs.
* `plot.TS()` provides MCMC diagnostic plots and summary plots.
* `plot.LDA_TS()` plots produce the combination of plots.

## Rodents data set
* Portal rodent data from [Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373) are now provided in a pre-formatted and ready-to-roll data object.
* Access the data using `data(rodents)`.
* Note, however, that the data in Christensen *et al.* 2018 are scaled according to trapping effort. The data included in LDATS are not, to allow for appropriate weighting. See [comparison vignette](https://weecology.github.io/LDATS/articles/paper-comparison.html) for further details.


## Comparison with [Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373)
* The [comparison vignette](https://weecology.github.io/LDATS/articles/paper-comparison.html) provides a step-by-step comparison of the LDATS pipeline to the analysis in Christensen *et al.* 2018. 
* The key differences are as follows:

      * The `document_term_table` in Christensen *et al.* 2018 was adjusted to account for variable trapping effort. The data included in LDATS are not adjusted, so that sampling periods can be weighted appropriately.
      * The LDA model selection criterion has changed (see LDA model AIC calculation, above), so that LDATS now identifies 6 topics compared to the 4 topics found in the paper.
      * LDATS will by default weight sampling periods according to the number of terms (see Document weighting, above). 
      * Despite these changes, the updated LDATS pipeline gives qualitatively similar results to the analysis in Christensen *et al.* 2018. 

# [LDATS 0.0.1](https://github.com/weecology/LDATS/commit/326506b9d7fb3e0223948d0245381963f83a2b37) 
*2017-11-16*

* Beginning initial development of package from [original code](https://github.com/emchristensen/Extreme-events-LDA) used in [Christensen *et al.* (2018)](https://doi.org/10.1002/ecy.2373).
