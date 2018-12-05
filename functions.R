library('rstan')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

envColors<-read.csv('envColor.csv',stringsAsFactors=FALSE)
envColors$col<-sprintf('#%s',envColors$color)
envCols<-structure(envColors$col,.Names=envColors$env)

stanCode<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
    //matrix<lower=0,upper=1>[nEnv,nRep] prop1;
    //matrix<lower=0,upper=1>[nEnv,nRep] prop2;
    matrix[nEnv,nRep] prop1;
    matrix[nEnv,nRep] prop2;
  }
  parameters {
    real metaFoldChange;
    real<lower=0> metaSd;
    vector[nRep] prop1BaseRaw[nEnv];
    vector[nEnv] metaRaw;
    real<lower=0> sigma;
    vector[nRep-1] repEffect;
    real<lower=0> baseSigma;
    real baseMu[nEnv];
  }
  transformed parameters{
    vector[nEnv] foldChange;
    vector[nRep] repEffectPlusZero;
    vector[nRep] prop1Base[nEnv];
    foldChange=metaFoldChange+metaRaw*metaSd;
    repEffectPlusZero[1]=0;
    repEffectPlusZero[2:nRep]=repEffect[1:(nRep-1)];
    for(ii in 1:nEnv){
      prop1Base[ii]=baseMu[ii]+repEffectPlusZero+prop1BaseRaw[ii]*baseSigma;
    }
  }
  model {
    for(ii in 1:nEnv){
      prop1BaseRaw[ii]~normal(0,1);
      prop1[ii,]~normal(prop1Base[ii],sigma);
      prop2[ii,]~normal(prop1Base[ii]+foldChange[ii],sigma);
    }
    //back~normal(-1.8,1);
    metaRaw~normal(0,1);
  }
'
mod <- stan_model(model_code = stanCode)

stanCode2<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
    matrix[nEnv,nRep] prop1;
    matrix[nEnv,nRep] prop2;
  }
  parameters {
    real metaFoldChange;
    real<lower=0> metaSd;
    real prop1Base[nEnv];
    vector[nEnv] metaRaw;
    real<lower=0> sigma;
    vector[nRep-1] repEffect1;
    vector[nRep-1] repEffect2;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    vector[nRep] repEffectPlusZero1;
    vector[nRep] repEffectPlusZero2;
    foldChange=metaFoldChange+metaRaw*metaSd;
    repEffectPlusZero1[1]=0;
    repEffectPlusZero1[2:nRep]=repEffect1;
    repEffectPlusZero2[1]=0;
    repEffectPlusZero2[2:nRep]=repEffect2;
  }
  model {
    for(ii in 1:nEnv){
      prop1[ii,]~normal(prop1Base[ii]+repEffectPlusZero1,sigma);
      prop2[ii,]~normal(prop1Base[ii]+repEffectPlusZero2+foldChange[ii],sigma);
    }
    //back~normal(-1.8,1);
    metaRaw~normal(0,1);
  }
'
mod2 <- stan_model(model_code = stanCode2)

plotRaw<-function(raw1,raw2,lab1=sub(' .*$','',colnames(raw1)[1]),lab2=sub(' .*$','',colnames(raw2)[1]),ylim=range(cbind(raw1,raw2)),cols=NULL){
  par(mfrow=c(2,ceiling(nrow(raw1)/2))) 
  notInCols<-rownames(raw1)[!rownames(raw1) %in% names(cols)]
  cols<-c(cols,structure(rep('black',length(notInCols)),.Names=notInCols))
  for(ii in 1:nrow(raw1)){
    plot(1,1,type='n',xlim=c(.5,2.5),ylim=ylim,log='y',yaxt='n',xlab='',ylab='Percent infected',main=rownames(raw1)[ii],xaxt='n')
    axis(1,1:2,c(lab1,lab2))
    segments(rep(1,ncol(raw1)),unlist(raw1[ii,]),rep(2,ncol(raw2)),unlist(raw2[ii,]))
    points(rep(1,ncol(raw1)),raw1[ii,],pch=21,bg=cols[rownames(raw1)[ii]])
    points(rep(2,ncol(raw2)),raw2[ii,],pch=21,bg=cols[rownames(raw1)[ii]])
    text(rep(1.5,ncol(raw1)),exp((log(raw1[ii,])+log(raw2[ii,]))/2),1:ncol(raw1))
    dnar::logAxis(las=1)
  }
}

findLims<-function(fit){
  mat<-as.matrix(fit)
  stats<-apply(mat,2,function(xx)c('lower'=unname(quantile(xx,c(.025))),'upper'=unname(quantile(xx,c(.975)))))
  selects<-stats[,grep('foldChange|metaFoldChange\\[',colnames(stats))]
  return(range(selects))
}

plotFit<-function(fit,virus,main='',cols=NULL,xlims=c(min(c(1,selects)),max(c(1,selects)))){
  mat<-as.matrix(fit)
  stats<-apply(mat,2,function(xx)c('mean'=mean(xx),'lower'=unname(quantile(xx,c(.025))),'upper'=unname(quantile(xx,c(.975)))))
  #selects<-exp(stats[,grep('foldChange|metaFoldChange|repEffect\\[',colnames(stats))])
  selects<-exp(stats[,grep('foldChange|metaFoldChange',colnames(stats))])
  foldIds<-grep('foldChange',colnames(selects))
  nums<-as.numeric(sub('.*\\[([0-9]+)\\]','\\1',colnames(selects)[foldIds]))
  colnames(selects)[foldIds]<-virus[nums]
  colnames(selects)[colnames(selects)=='metaFoldChange']<-'Overall'
  colnames(selects)[colnames(selects)=='repEffect[1]']<-'Replicate 2'
  colnames(selects)[colnames(selects)=='repEffect[2]']<-'Replicate 3'
  selects<-selects[,ncol(selects):1]
  notInCols<-colnames(selects)[!colnames(selects) %in% names(cols)]
  cols<-c(cols,structure(rep('black',length(notInCols)),.Names=notInCols))
  plot(1,1,ylim=c(1,ncol(selects)),xlim=xlims,log='x',type='n',ylab='',xlab='Fold change',main=main,yaxt='n',xaxt='n',mgp=c(2.5,.8,0))
  axis(2,1:ncol(selects),colnames(selects),las=1)
  segments(selects['lower',],1:ncol(selects),selects['upper',],1:ncol(selects))
  points(selects['mean',],1:ncol(selects),pch=21,bg=cols[colnames(selects)])
  dnar::logAxis(1)
  if(xlims[2]<3&xlims[1]>.1)axis(1,.5)
  if(xlims[2]<10&xlims[1]>.5)axis(1,.7)
  abline(v=1,lty=3)
}

compareAlleles<-function(mod,allele1,allele2,dat1,dat2=dat1){
  message(allele1,' - ',allele2)
  dat=list(
    nEnv=nrow(dat1),
    nRep=3,
    prop1=log(dat1[,grep(sprintf('^%s [0-9]+$',allele1),colnames(dat1))]/100),
    prop2=log(dat2[,grep(sprintf('^%s [0-9]+$',allele2),colnames(dat2))]/100)
  )
  fit <- sampling(mod, data = dat, iter=50000, chains=30,thin=25,control=list(adapt_delta=.9))
}



