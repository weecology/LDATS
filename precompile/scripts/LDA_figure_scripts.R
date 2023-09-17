# Scripts for making figures from LDA analyses

library(ggplot2)
library(gridExtra)
library(dplyr)


#cbPalette <- c( "#E69F00","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")
cbPalette <- c( "#e19c02","#999999", "#56B4E9", "#0072B2", "#D55E00", "#F0E442", "#009E73", "#CC79A7")




#' Plot gamma
#' 
#' This plotting function plots the component communities over time
#' It's used by the plot_component_communities function, but can also take
#' any data frame as input as long as it's the right form -- used by the simulations
#' 
#' @param gamma_frame a data frame containing columns for date, relabund, and community
#' @param ntopics number of topics
#' @param ylab label for y axis (optional)
#' 
#' @return a ggplot object

ltypes = c('solid','longdash','longdash','solid','solid')


plot_gamma = function(gamma_frame,ntopics,ylab='',colors=cbPalette) {
  g = ggplot(gamma_frame, aes(x=date,y=relabund,colour=community)) + 
    #geom_point() +
    geom_line(aes(size=community)) +
    scale_size_manual(values=c(1,1,1,1,1), guide=FALSE) +
    scale_y_continuous(name=ylab,limits=c(0,1)) +
    scale_x_date(name='') +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          #panel.background = element_blank(),
          panel.border=element_rect(colour='black',fill=NA),
          legend.position='none') +
    scale_colour_manual(name="Component\nCommunity",
                        breaks=as.character(seq(ntopics)),
                        values=colors[1:ntopics],
                        guide=FALSE)
  return(g)
}


#' Plot component communities
#' 
#' Plots timeseries of component communities (topics)
#' 
#' @param ldamodel object of class LDA_VEM created by the function LDA in topicmodels package
#' @param ntopics number of topics used in ldamodel
#' @param xticks vector of dates for x-axis labels
#' @param ylab y axis label (optional)
#' @param topic_order order of topics (for color control)
#' 
#' @return ggplot object
#' 
#' @example plot_component_communities(ldamodel,ntopics,period_dates$date)

plot_component_communities = function(ldamodel,ntopics,xticks,ylab='',topic_order = seq(ntopics),colors = cbPalette) {
  
  z = posterior(ldamodel)
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (t in topic_order) {
    ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=z$topics[,t],community = as.factor(rep(t,length(z$topics[,1])))))
  }
  g = plot_gamma(ldaplot,ntopics,ylab,colors)
  return(g) 
}


#' Plot histogram of changepoint locations -- WIP
#' 
#' @param results results object from changepoint_model function
#' @param year_continuous vector of dates/xaxis units
#' 
#' @return ggplot object

chpoint_histogram = function(results,year_continuous) {
  npts = dim(results$saved)[1]
  nrep = dim(results$saved)[3]
  df = as.data.frame(t(results$saved[,1,]))

  ggplot(data = df, aes(x=value)) +
    for (n in seq(npts)) {
      geom_histogram(data=data.frame(value=results$saved[n,1,]))
    }
    geom_histogram(aes(y=..count../sum(..count..),col=color),
                   binwidth = .25) +
    labs(x='') +
    xlim(range(year_continuous))
    theme(axis.text=element_text(size=12),
          panel.border=element_rect(colour='black',fill=NA))
  return(h)
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
#' @param topic_order order of topics -- for making this bar graph relate to the component community graph
#' 
#' @return barplots of the n component communities
#'
#' @example plot_community_composition(community_composition(ldamodel))

plot_community_composition = function(composition,topic_order=1:dim(composition)[1],colors = cbPalette) {
  nspecies = dim(composition)[2]
  topics = dim(composition)[1]
  ylimits = c(0,round(max(composition),1)+.1)
  par(mfrow=c(1,topics))
  j=1
  for (i in topic_order) {
    x = barplot(composition[i,],ylim=ylimits,col=colors[i],main=paste('Community-type',j),las=2)
    j=j+1
  }
  par(mfrow=c(1,1))
}

#' ggplot version of plot_community_composition
#' 
#' @param composition matrix of species composition of topics; as in output of community_composition()
#' @param topic_order order of topics -- for making this bar graph relate to the component community graph
#' @param ylim vector of limits for yaxis
#' @param colors color palette specification
#' @param title T/F if you want title or not
#' @param ylabels T/F if you want yaxis label or not
#' 
#' @return barplots of the n component communities
#' 
#' 
#' 
plot_community_composition_gg = function(composition,topic_order,ylim,colors=cbPalette,title=T,ylabels = T) {
  topics = dim(composition)[1]
  community = c()
  for (j in 1:topics) {community=append(community,rep(j,length(composition[j,])))}
  relabund = c()
  for (j in 1:topics) {relabund=append(relabund,composition[j,])}
  species=c()
  for (j in 1:topics) {species=append(species,colnames(composition))}
  comp = data.frame(community = community,relabund=relabund,species=factor(species, levels = colnames(composition)))
  grass = filter(comp,species %in% c('BA','PH','DO','DS','PF','PL','PM','RF','RO','RM','SH','SF','SO'))
  p = list()
  j = 1
  for (i in topic_order) {
    if (j == 1 && ylabels == T) {ylabel='% Composition'} else {ylabel=''}
    x <- ggplot(data=comp[comp$community==i,], aes(x=species, y=relabund)) +
      geom_bar(stat='identity',fill=colors[i])  +
      geom_bar(data=grass[grass$community==i,],aes(x=species,y=relabund),fill=colors[i],stat='identity',alpha=0,size=.5,color='black') +
        theme(axis.text=element_text(size=9),
              panel.background = element_blank(),
              panel.border=element_rect(colour='black',fill=NA),
              axis.text.x = element_text(angle = 90,hjust=0,vjust=.5,size=6),
              plot.margin = unit(c(0,1,0,0),"mm"),
              axis.text.y = element_text(angle=0,size=8,vjust=.5,hjust=.5),
              plot.title = element_text(hjust = 0.5,size=10.5)) +
      scale_x_discrete(name='') +
      scale_y_continuous(name=ylabel,limits = ylim) +
      geom_hline(yintercept = 0)  +
      if (title==T) {ggtitle(paste('Community-\ntype',j))} else {ggtitle('')}

    p[[j]] <- x
    j=j+1
  }
 
  return(p)
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
  
  theta1=matrix(apply(results$theta,2,mean),ntopics,length(xticks))
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (t in seq(ntopics)) {
    ldaplot = rbind(ldaplot,data.frame(date=xticks,relabund=theta1[t,],community = as.factor(rep(t,length(xticks)))))
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
    scale_colour_manual(name="Component \nCommunity",
                        breaks=as.character(seq(ntopics)),
                        values=cbPalette[1:ntopics])
}




#' Plot component communities including credible intervals --- WIP!
#' so far only works with 3 topics -- need to make more general
#' 
#' Plots timeseries of component communities (topics) from LDA using Gibbs sampler
#' 
#' @param ldamodel model output of LDA model using Gibbs
#' @param ntopics number of topics used in ldamodel
#' @param xticks vector of dates for x-axis labels
#' 
#' @return None
#' 
#' @example plot_component_communities_gibbs_credible(ldamodel,ntopics,period_dates$date)


plot_component_communities_gibbs_credible = function(ldamodel,ntopics,xticks) {
  nsteps = length(xticks)
  theta1=matrix(apply(ldamodel$theta,2,mean),ntopics,nsteps)
  thetasd = matrix(apply(ldamodel$theta,2,sd),ntopics,nsteps)
  
  thetadf = data.frame(grp1 = theta1[1,],grp2 = theta1[2,],sd1=thetasd[1,],sd2=thetasd[2,],grp3=theta1[3,],sd3=thetasd[3,])
  ggplot(thetadf) +
    geom_ribbon(data = thetadf, mapping = aes_string(x = xticks, ymin = thetadf$grp1-1.96*thetadf$sd1, ymax = thetadf$grp1+1.96*thetadf$sd1,alpha=.4), fill = "grey") +
    geom_line(aes(y = thetadf$grp1,x=xticks)) +
    geom_ribbon(data = thetadf, mapping = aes_string(x = xticks, ymin = thetadf$grp2-1.96*thetadf$sd2, ymax = thetadf$grp2+1.96*thetadf$sd2,alpha=.4), fill = "pink") +
    geom_line(aes(y = thetadf$grp2,x=xticks),color='red') +
    geom_ribbon(data = thetadf, mapping = aes_string(x = xticks, ymin = thetadf$grp3-1.96*thetadf$sd3, ymax = thetadf$grp3+1.96*thetadf$sd3,alpha=.4), fill = "lightgreen") +
    geom_line(aes(y = thetadf$grp3,x=xticks),color='forestgreen') +
    scale_y_continuous(name='Relative Abundance',limits=c(0,1)) +
    scale_x_date(name='')


}
