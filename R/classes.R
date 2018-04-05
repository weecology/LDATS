setMethod("plot", signature(x = "LDA", y = "ANY"), 
  function(x, y, cols) {
    plot_LDA(x, cols)
  }
)