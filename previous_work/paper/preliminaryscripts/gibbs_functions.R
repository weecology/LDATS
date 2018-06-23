library(Rcpp)

gibbs.samp=function(dat.agg,ngibbs,ncommun,a.betas,a.theta){
  
  arv=apply(dat.agg,1,sum)
  nplots=nrow(dat.agg)
  ntrees=max(arv)
  nspp=ncol(dat.agg)
  
  dat=matrix(0,nplots,ntrees)
  for (i in 1:nplots){
    z=rep(1:nspp,dat.agg[i,])
    dat[i,1:length(z)]=z
  }
  
  ind.na=which(dat==0)
  ind.notna=which(dat!=0)
  dat[ind.na]=-999
  
  #--------------------------------------
  arv1=matrix(arv,ncommun,nplots,byrow=T)
  nspp=max(dat,na.rm=T)
  aux.plot=matrix(1:nplots,nplots,ntrees)
  
  z=matrix(NA,nplots,ntrees) #create indicator of which community each tree belongs to
  z[]=sample(1:ncommun,size=nplots*ntrees,replace=T)
  z[ind.na]=-999
  
  nji=matrix(0,ncommun,nspp)
  nj=matrix(NA,ncommun,1)
  for (i in 1:ncommun){
    ind=which(z==i)
    nj[i]=length(ind)
    dat1=dat[ind]
    dat2=data.frame(table(dat1))
    nji[i,as.numeric(as.character(dat2$dat1))]=dat2$Freq
  }
  
  ndj=matrix(0,ncommun,nplots)
  for (i in 1:nplots){
    z1=data.frame(table(z[i,]))
    z1$Var1=as.numeric(as.character(z1$Var1))
    z2=z1[z1$Var1!=-999,]
    ndj[z2$Var1,i]=z2$Freq
  }
  
  #gibbs
  param=list(z=z,ndj=ndj,nji=nji,nj=nj)
  seq1=floor(seq(from=100,to=ngibbs,length.out=100))
  store.theta=matrix(NA,100,length(ndj))
  store.beta=matrix(NA,100,length(nji))
  
  param$dat=dat
  param$a.betas=a.betas
  param$a.theta=a.theta
  param$nspp=nspp
  param$ntrees=ntrees
  
  thin=3
  logL=-Inf
  logL.new=rep(NA,ngibbs)
  oo=1
  for (i in 1:ngibbs){
    print(i)
    k=lda(ind.notna, param$z, param$nj, param$dat, param$nji, param$ndj, aux.plot, param$a.betas, param$a.theta, nspp, ncommun,thin)
    param$z=k$z
    param$nj=k$nj
    param$nji=k$nji
    param$ndj=k$ndj
    
    #create theta and beta
    k=param$nji+1
    gammas=matrix(rgamma(ncommun*nspp,k,1),ncommun,nspp) #dirichlet generation
    soma=matrix(gammas%*%rep(1,nspp),ncommun,nspp)
    beta1=gammas/soma
    
    k=param$ndj+1
    gammas=matrix(rgamma(ncommun*nplots,k,1),ncommun,nplots) #dirichlet generation
    soma=matrix(rep(1,ncommun)%*%gammas,ncommun,nplots,byrow=T)
    theta1=gammas/soma
    
    if (i%in%seq1){
      store.theta[oo,]=theta1 
      store.beta[oo,]=beta1
      oo=oo+1
    }
    
    #loglikelihood
    exp1=t(theta1)%*%beta1
    exp2=dat.agg*log(exp1) #logL of multinomial distribution
    logL.new[i]=sum(exp2)
  }
  
  list('logL'=logL.new,'beta'=store.beta,'theta'=store.theta)
}

cppFunction('
  Rcpp::List lda(IntegerVector seq1, IntegerMatrix z, IntegerVector nj, IntegerMatrix dat, IntegerMatrix nji, IntegerMatrix ndj, IntegerMatrix auxplot, double abetas, double atheta, int nspp, int ncommun, int thin) {
    int n = seq1.size();
    IntegerMatrix z1 = z;
    IntegerVector nj1 = nj;
    IntegerMatrix nji1 = nji;
    IntegerMatrix ndj1 = ndj; 

    for (int dd = 0; dd < thin; dd++){
      for (int i = 0; i < n; i++) {
        int ind = seq1[i]-1; //vectors start at 0
        if (z1[ind] > 0){
          int sp = dat[ind]-1;
          int comm = z1[ind]-1;
          int plot = auxplot[ind]-1;
        
          nj1[comm] = nj1[comm]-1;
          nji1(comm,sp) = nji1(comm,sp)-1;
          ndj1(comm,plot) = ndj1(comm,plot)-1;
            
          NumericVector prob=log(nji1(_,sp)+abetas) + log(ndj1(_,plot)+atheta) - log(nj1+nspp*abetas);
          NumericVector prob1=prob-max(prob);
          NumericVector prob2=exp(prob1);
          NumericVector prob3=prob2/sum(prob2);
          
          //multinomial distribution
          double res = 0;
          double teste = runif(1,0,1)[0];
          int fim = 0;
          for (int j = 0; j < ncommun; j++){
            res = res + prob3[j];
            if (teste < res) {
              fim=j;
              break;
            }
          }
            
          z1[ind] = fim+1;
          nj1[fim] = nj1[fim]+1;
          nji1(fim,sp) = nji1(fim,sp)+1;
          ndj1(fim,plot) = ndj1(fim,plot)+1;
        }
      }
    }
    return (Rcpp::List::create(Rcpp::Named("z") = z1, Rcpp::Named("nj") = nj1, Rcpp::Named("nji") = nji1, Rcpp::Named("ndj") = ndj1));
  }
')