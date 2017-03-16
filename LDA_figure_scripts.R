# Scripts for making figures from LDA analyses

plot_component_communities = function(ldamodel,k,xticks) {
  # function to take output of LDA function and plot time series of component communities
  # Inputs:
  #   ldamodel = output of function LDA
  #   k = number of topics
  #   xticks = index of x-axis, e.g. a vector of dates
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  z = posterior(ldamodel)
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (k1 in seq(k)) {
    ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=z$topics[,k1],community = as.factor(rep(k1,length(z$topics[,1])))))
  }
  ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_point() +
    geom_line(size=1) +
    scale_y_continuous(name='Percent Similarity') +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    scale_colour_manual(name="",
                        breaks=as.character(seq(k)),
                        values=cbPalette[1:k])
  
}

plot_component_communities_smooth = function(ldamodel,k,xticks) {
  # function to take output of LDA function and plot time series of component communities
  # Inputs:
  #   ldamodel = output of function LDA
  #   k = number of topics
  #   xticks = index of x-axis, e.g. a vector of dates
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  z = posterior(ldamodel)
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (k1 in seq(k)) {
    ldacomm = data.frame(date=xticks,relabund=z$topics[,k1],community = as.factor(rep(k1,length(z$topics[,1]))))
    ldacomm$smooth = rep(NA)
    for (n in seq(length(ldacomm$smooth)-5)) {
      ldacomm$smooth[n+2] = mean(ldacomm$relabund[n:(n+5)])
    }
    ldaplot = rbind(ldaplot,ldacomm)
  }
  
  ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_point() +
    geom_line(size=1,aes(y=smooth)) +
    scale_y_continuous(name='Percent Similarity') +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    scale_colour_manual(name="",
                        breaks=as.character(seq(k)),
                        values=cbPalette[1:k])
  
}

community_composition = function(ldamodel) {
  # function that returns a matrix of the species composition of the component communities (to 3 decimal places)
  return(structure(round(exp(ldamodel@beta), 3), dimnames = list(NULL, ldamodel@terms)))
}

