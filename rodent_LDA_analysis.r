# LDA analysis (with figures) for the extreme events/ LDA project

# LDA =================================
library(topicmodels)
library(ggplot2)
set.seed(1)

# load csv of abundance of each sp per plot
dat = read.csv('Rodent_table_dat.csv',na.strings = '',as.is=T)

# dates of trapping periods
perdat= read.csv('period_dates_single.csv')

perdat$date = as.Date(perdat$date,format='%m/%d/%Y')


d = ymd(as.Date(perdat$date[1:length(dat[,1])]))
year_continuous = 1970 + as.integer(julian(d)) / 365.25
x = data.frame(
  year_continuous = year_continuous,
  sin_year = sin(year_continuous * 2 * pi),
  cos_year = cos(year_continuous * 2 * pi)
)



# LDA models: groups from 2 to 5
nstart = 20 # For the final analysis, maybe do 1000
ldamodel2 = LDA(dat,2,control=list(estimate.alpha=F,alpha=1, nstart = nstart),method="VEM")
ldamodel3 = LDA(dat,3,control=list(estimate.alpha=F,alpha=1, nstart = nstart),method="VEM")
ldamodel4 = LDA(dat,4,control=list(estimate.alpha=F,alpha=1, nstart = nstart),method="VEM")
ldamodel5 = LDA(dat,5,control=list(estimate.alpha=F,alpha=1, nstart = nstart),method="VEM")

# composition of component communities, to identify most important species
structure(round(exp(ldamodel4@beta), 3), dimnames = list(NULL, ldamodel2@terms))


# ===========================================
# figures

dates = as.Date(perdat$date[1:length(dat[,1])])

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_component_communities = function(ldamodel,k) {
  # function to take output of LDA function and plot time series of component communities
  z = posterior(ldamodel)
  ldaplot = data.frame(date=c(),relabund = c(), community = c())
  for (k1 in seq(k)) {
    ldaplot = rbind(ldaplot,data.frame(date=dates,relabund=z$topics[,k1],community = as.factor(rep(k1,436))))
  }
  ggplot(ldaplot, aes(x=date,y=relabund,colour=community)) + 
    geom_point() +
    geom_line(size=1) +
    scale_x_date(name='') +
    scale_y_continuous(name='Percent Similarity') +
    theme(axis.text=element_text(size=18),
          axis.title=element_text(size=18),
          legend.text=element_text(size=18),
          legend.title=element_text(size=18)) +
    scale_colour_manual(name="",
                        breaks=as.character(seq(k)),
                        values=cbPalette[1:k])

}


plot_component_communities(ldamodel2,2)
plot_component_communities(ldamodel3,3)
plot_component_communities(ldamodel4,4)
plot_component_communities(ldamodel5,5)
