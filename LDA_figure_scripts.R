# Scripts for making figures from LDA analyses

library(ggplot2)







#' Plot component communities
#' 
#' Plots timeseries of component communities (topics)
#' 
#' @param ldamodel object of class LDA_VEM created by the function LDA in topicmodels package
#' @param ntopics number of topics used in ldamodel
#' @param xticks vector of dates for x-axis labels
#' 
#' @return None
#' 
#' @example plot_component_communities(ldamodel,ntopics,period_dates$date)

plot_component_communities = function(ldamodel,ntopics,xticks) {
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  z = posterior(ldamodel)
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (t in seq(ntopics)) {
    ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=z$topics[,t],community = as.factor(rep(t,length(z$topics[,1])))))
  }
  ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_point() +
    geom_line(size=1) +
    scale_y_continuous(name='Relative Abundance',limits=c(0,1)) +
    scale_x_date(name='') +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    scale_colour_manual(name="Component\nCommunity",
                        breaks=as.character(seq(ntopics)),
                        values=cbPalette[1:ntopics])
  
}



#' Plot component communities -- smoothed
#' 
#' Plots timeseries of component communities (topics)
#' Smooths using a simple moving window average
#' Plots raw data as dots, smoothed as lines
#' 
#' @param ldamodel object of class LDA_VEM created by the function LDA in topicmodels package
#' @param ntopics number of topics used in ldamodel
#' @param xticks vector of dates for x-axis labels
#' @param smooth_factor size of moving window average -- higher value is smoother
#' 
#' @return None
#' 
#' @example plot_component_communities_smooth(ldamodel,ntopics,period_dates$date,5)

plot_component_communities_smooth = function(ldamodel,ntopics,xticks,smooth_factor) {
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  z = posterior(ldamodel)
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (t in seq(ntopics)) {
    ldacomm = data.frame(date=xticks,relabund=z$topics[,t],community = as.factor(rep(t,length(z$topics[,1]))))
    ldacomm$smooth = rep(NA)
    for (n in seq(length(ldacomm$smooth)-smooth_factor)) {
      ldacomm$smooth[n+floor(smooth_factor/2)] = mean(ldacomm$relabund[n:(n+smooth_factor)])
    }
    ldaplot = rbind(ldaplot,ldacomm)
  }
  
  ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_point() +
    geom_line(size=1,aes(y=smooth)) +
    scale_y_continuous(name='Relative Abundance') +
    scale_x_date(name='') +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    scale_colour_manual(name="Component \nCommunity",
                        breaks=as.character(seq(ntopics)),
                        values=cbPalette[1:ntopics])
  
}


#' Make table of species composition of topics
#' 
#' @param ldamodel  object of class LDA_VEM created by the function LDA in topicmodels package
#' 
#' @return table of species composition of the topics in ldamodel, 3 decimal places
#' 
#' @example community_composition(ldamodel)

community_composition = function(ldamodel) {
  return(structure(round(exp(ldamodel@beta), 3), dimnames = list(NULL, ldamodel@terms)))
}


#' Plot species composition of topics
#' 
#' @param composition matrix of species composition of topics; as in output of community_composition()
#' @param ylimits vector of (ymin,ymax) for plotting species composition
#' 
#' @return None
#'
#' @example plot_community_composition(community_composition(ldamodel),c(0,1))

plot_community_composition = function(composition,ylimits) {
  nspecies = dim(composition)[2]
  topics = dim(composition)[1]
  par(mfrow=c(topics,1))
  for (i in 1:topics) {plot(1:nspecies,composition[i,],type='l',xlab='',ylim=ylimits,col=i,lwd=2,xaxt='n',ylab='')
    axis(1,at=1:nspecies,labels=colnames(composition))}
}


#' Plot component communities
#' 
#' Plots timeseries of component communities (topics) from LDA using Gibbs sampler
#' 
#' @param ldamodel model output of LDA model using Gibbs
#' @param ntopics number of topics used in ldamodel
#' @param xticks vector of dates for x-axis labels
#' 
#' @return None
#' 
#' @example plot_component_communities_gibbs(ldamodel,ntopics,period_dates$date)

plot_component_communities_gibbs = function(results,ntopics,xticks) {
  
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  theta1=matrix(apply(results$theta,2,mean),ntopics,length(xticks))
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (t in seq(ntopics)) {
    ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=theta1[t,],community = as.factor(rep(t,length(xticks)))))
  }
  
  ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_point() +
    geom_line(size=1) +
    scale_y_continuous(name='Percent Similarity',limits=c(0,1)) +
    scale_x_date(name='') +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    scale_colour_manual(name="",
                        breaks=as.character(seq(ntopics)),
                        values=cbPalette[1:ntopics])
}
