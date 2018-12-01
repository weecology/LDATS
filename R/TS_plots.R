
plot.TS_fit <- function(x, ..., plot_type = "diagnostic"){
  if (plot_type == "diagnostic"){
    TS_diagnostics_plot(x, ...)
  } else if (plot_type == "summary"){
    TS_summary_plot(x, ...)
  }
}


TS_summary_plot <- function(x, ...){

}

rho_hist <- function(x){

}