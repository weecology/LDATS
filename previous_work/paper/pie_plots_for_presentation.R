

dat = read.csv('Rodent_table_dat.csv',as.is=T)

# ========================================================================================
# pie plot of species composition for a particular period


lbls = c('B. taylori','D. merriami','D. ordii','D. spectabilis','N. albigula','O. leucogaster','O. torridus',
         'C. baileyi','P. eremicus','P. flavus','C. hispidus','C. intermedius','P. leucopus','P. maniculatus',
         'C. penicillatus','R. fulvescens','R. montanus','R. megalotis','S. fulviventer','S. hispidus','S. ochrognathus')

ppcols = c('#F79646','#4BACC6','#8064A2','#9BBB59','#C0504D','#4F81BD','#1F497D','#EEECE1','#D9CC25','#FCD5B5',
           '#FCD5B5','#B7DEE8','#CCC1DA','#D7E4BD','#E6B9B8','#B9CDE5','#8EB4E3','#C4BD97','#7F7F7F','#000000',
           '#1F497D')

# period  for presentation, used 426 (May 2015) and 22 (April 1979)
n = 22


pie(unlist(dat[n,]),
    labels=NA,
    clockwise=TRUE,
    col=ppcols,
    border="white",
    radius=0.7,
    cex=0.8)
legend("right",legend=lbls[dat[n,]>0],bty="n",
       fill=ppcols[dat[n,]>0])
       