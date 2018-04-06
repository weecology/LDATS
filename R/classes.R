setMethod("plot", signature(x = "LDA", y = "ANY"), 
  function(x, y, cols = NULL) {
    plot_LDA(x, cols = cols)
  }
)