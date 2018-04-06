setMethod("plot", signature(x = "LDA", y = "ANY"), 
  function(x, y) {
    plot_LDA(x, cols = NULL)
  }
)