library('rstan')
if(!dir.exists('out'))dir.create('out')
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

stanCode3<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
    matrix[nEnv,nRep] prop1;
    matrix[nEnv,nRep] prop2;
    real wtSd;
    real wtMu;
  }
  parameters {
    real metaFoldChange;
    real<lower=0> metaSd;
    real prop1Base[nEnv];
    vector[nEnv-1] metaRaw;
    real<lower=0> sigma;
    real wtFoldChange;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    foldChange[1]=wtFoldChange;
    foldChange[2:nEnv]=metaFoldChange+metaRaw*metaSd;
  }
  model {
    prop1[1,]~normal(prop1Base[1],sigma);
    prop2[1,]~normal(prop1Base[1]+foldChange[1],sigma);
    for(ii in 2:nEnv){
      prop1[ii,]~normal(prop1Base[ii],sigma);
      prop2[ii,]~normal(prop1Base[ii]+foldChange[1]+foldChange[ii],sigma);
    }
    metaRaw~normal(0,1);
    wtFoldChange~normal(wtMu,wtSd);
  }
'
mod3 <- stan_model(model_code = stanCode3)

stanCode4<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
    matrix<lower=0,upper=1>[nEnv,nRep] prop1;
    matrix<lower=0,upper=1>[nEnv,nRep] prop2;
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
    matrix<upper=0>[nEnv,nRep] back1;
    matrix<upper=0>[nEnv,nRep] back2;
    matrix[nEnv,nRep] trueProp1Raw;
    matrix[nEnv,nRep] trueProp2Raw;
    matrix<lower=2000>[nEnv,nRep] nCells1;
    matrix<lower=2000>[nEnv,nRep] nCells2;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    vector[nRep] repEffectPlusZero;
    vector[nRep] prop1Base[nEnv];
    matrix<lower=0,upper=1>[nEnv,nRep] propPlusBack1;
    matrix<lower=0,upper=1>[nEnv,nRep] propPlusBack2;
    foldChange=metaFoldChange+metaRaw*metaSd;
    repEffectPlusZero[1]=0;
    repEffectPlusZero[2:nRep]=repEffect[1:(nRep-1)];
    for(ii in 1:nEnv){
      prop1Base[ii]=baseMu[ii]+repEffectPlusZero+prop1BaseRaw[ii]*baseSigma;
      propPlusBack1[ii,]=inv_logit(trueProp1Raw[ii,]*sigma+to_row_vector(prop1Base[ii]))+inv_logit(back1[ii,]);
      propPlusBack2[ii,]=inv_logit(trueProp2Raw[ii,]*sigma+to_row_vector(prop1Base[ii])+foldChange[ii])+inv_logit(back2[ii,]);
      for(jj in 1:nRep){
        propPlusBack1[ii,jj]=fmin(propPlusBack1[ii,jj],.9999);
        propPlusBack2[ii,jj]=fmin(propPlusBack2[ii,jj],.9999);
      }
    }
  }
  model {
    baseSigma~gamma(1,.1);
    sigma~gamma(1,.1);
    metaRaw~normal(0,1);
    for(ii in 1:nEnv){
      prop1BaseRaw[ii]~normal(0,1);
      //trueProp1[ii,]~normal(prop1Base[ii],sigma);
      //trueProp2[ii,]~normal(prop1Base[ii]+foldChange[ii],sigma);
      trueProp1Raw[ii,]~normal(0,1);
      trueProp2Raw[ii,]~normal(0,1);
      back1[ii,]~normal(-6.43,.97);
      back2[ii,]~normal(-6.43,.97);
      nCells1[ii,]~normal(5000,750);
      nCells2[ii,]~normal(5000,750);
      prop1[ii,]~normal(propPlusBack1[ii,],sqrt(propPlusBack1[ii,].*(1-propPlusBack1[ii,])./nCells1[ii,]));
      prop2[ii,]~normal(propPlusBack2[ii,],sqrt(propPlusBack2[ii,].*(1-propPlusBack2[ii,])./nCells2[ii,]));
    }
  }
'
stanCode5<-'
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
    real baseMuRaw[nEnv];
    matrix<upper=0>[nEnv,nRep] back1;
    matrix<upper=0>[nEnv,nRep] back2;
    real metaBase;
    real<lower=0>metaBaseSd;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    vector[nRep] repEffectPlusZero;
    vector[nRep] prop1Base[nEnv];
    matrix<lower=0>[nEnv,nRep] multBack1;
    matrix<lower=0>[nEnv,nRep] multBack2;
    foldChange=metaFoldChange+metaRaw*metaSd;
    repEffectPlusZero[1]=0;
    repEffectPlusZero[2:nRep]=repEffect[1:(nRep-1)];
    for(ii in 1:nEnv){
      prop1Base[ii]=metaBase+metaBaseSd*baseMuRaw[ii]+repEffectPlusZero+prop1BaseRaw[ii]*baseSigma;
      multBack1[ii,]=log((exp(back1[ii,])+exp(to_row_vector(prop1Base[ii])))./exp(to_row_vector(prop1Base[ii])));
      multBack2[ii,]=log((exp(back2[ii,])+exp(to_row_vector(prop1Base[ii]+foldChange[ii])))./exp(to_row_vector(prop1Base[ii]+foldChange[ii])));
    }
  }
  model {
    for(ii in 1:nEnv){
      prop1BaseRaw[ii]~normal(0,1);
      prop1[ii,]~normal(prop1Base[ii]+to_vector(multBack1[ii,]),sigma);
      prop2[ii,]~normal(prop1Base[ii]+foldChange[ii]+to_vector(multBack2[ii,]),sigma);
      back1[ii,]~normal(-6.43,.97);
      back2[ii,]~normal(-6.43,.97);
    }
    baseMuRaw~normal(0,1);
    metaRaw~normal(0,1);
    sigma~gamma(1,.1);
    baseSigma~gamma(1,.1);
    metaSd~gamma(1,.1);
  }
'
mod4 <- stan_model(model_code = stanCode5)

stanCode6<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
    matrix[nEnv,nRep] prop1;
    matrix[nEnv,nRep] prop2;
  }
  parameters {
    real metaFoldChange;
    real<lower=0> metaSd;
    real prop1BaseRaw[nEnv];
    vector[nEnv] metaRaw;
    real<lower=0> sigma;
    vector[nRep-1] repEffect1;
    vector[nRep-1] repEffect2;
    matrix<upper=0>[nEnv,nRep] back1;
    matrix<upper=0>[nEnv,nRep] back2;
    real metaBase;
    real<lower=0>metaBaseSd;
    real baseMuRaw[nEnv];
    real<lower=0> baseSigma;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    vector[nRep] repEffectPlusZero1;
    vector[nRep] repEffectPlusZero2;
    vector[nRep] prop1Base[nEnv];
    vector[nRep] prop2Base[nEnv];
    matrix<lower=0>[nEnv,nRep] multBack1;
    matrix<lower=0>[nEnv,nRep] multBack2;
    foldChange=metaFoldChange+metaRaw*metaSd;
    repEffectPlusZero1[1]=0;
    repEffectPlusZero1[2:nRep]=repEffect1;
    repEffectPlusZero2[1]=0;
    repEffectPlusZero2[2:nRep]=repEffect2;
    for(ii in 1:nEnv){
      prop1Base[ii]=metaBase+metaBaseSd*baseMuRaw[ii]+repEffectPlusZero1+prop1BaseRaw[ii]*baseSigma;
      prop2Base[ii]=metaBase+metaBaseSd*baseMuRaw[ii]+repEffectPlusZero2+prop1BaseRaw[ii]*baseSigma+foldChange[ii];
      multBack1[ii,]=log((exp(back1[ii,])+exp(to_row_vector(prop1Base[ii])))./exp(to_row_vector(prop1Base[ii])));
      multBack2[ii,]=log((exp(back2[ii,])+exp(to_row_vector(prop2Base[ii])))./exp(to_row_vector(prop2Base[ii])));
    }
  }
  model {
    for(ii in 1:nEnv){
      prop1BaseRaw[ii]~normal(0,1);
      prop1[ii,]~normal(prop1Base[ii]+to_vector(multBack1[ii,]),sigma);
      prop2[ii,]~normal(prop2Base[ii]+to_vector(multBack2[ii,]),sigma);
      back1[ii,]~normal(-6.43,.97);
      back2[ii,]~normal(-6.43,.97);
    }
    baseMuRaw~normal(0,1);
    metaRaw~normal(0,1);
    sigma~gamma(1,.1);
    baseSigma~gamma(1,.1);
    metaSd~gamma(1,.1);
  }
'
mod5 <- stan_model(model_code = stanCode6)

stanCode7<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
    matrix[nEnv,nRep] prop1;
    matrix[nEnv,nRep] prop2;
    real wtSd;
    real wtMu;
  }
  parameters {
    real metaFoldChange;
    real<lower=0> metaSd;
    real prop1Base[nEnv];
    vector[nEnv-1] metaRaw;
    real<lower=0> sigma;
    real wtFoldChange;
    matrix<upper=0>[nEnv,nRep] back1;
    matrix<upper=0>[nEnv,nRep] back2;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    matrix<lower=0>[nEnv,nRep] multBack1;
    matrix<lower=0>[nEnv,nRep] multBack2;
    foldChange[1]=wtFoldChange;
    foldChange[2:nEnv]=metaFoldChange+metaRaw*metaSd;
    for(ii in 1:nEnv){
      multBack1[ii,]=log((exp(back1[ii,])+exp(prop1Base[ii]))./exp(prop1Base[ii]));
      multBack2[ii,]=log((exp(back2[ii,])+exp(prop1Base[ii]+foldChange[1]+foldChange[ii]))./exp(prop1Base[ii]+foldChange[1]+foldChange[ii]));
    }
  }
  model {
    prop1[1,]~normal(prop1Base[1]+multBack1[1,],sigma);
    prop2[1,]~normal(prop1Base[1]+foldChange[1]+multBack2[1,],sigma);
    for(ii in 2:nEnv){
      prop1[ii,]~normal(prop1Base[ii]+multBack1[ii,],sigma);
      prop2[ii,]~normal(prop1Base[ii]+foldChange[1]+foldChange[ii]+multBack2[ii,],sigma);
      back1[ii,]~normal(-6.43,.97);
      back2[ii,]~normal(-6.43,.97);
    }
    metaRaw~normal(0,1);
    wtFoldChange~normal(wtMu,wtSd);
  }
'
mod6 <- stan_model(model_code = stanCode7)





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
  stats<-pullRanges(fit)
  selects<-stats[,grep('foldChange|metaFoldChange\\[',colnames(stats))]
  return(range(selects))
}
pullRanges<-function(fit,convertFunc=function(yy)yy){
  mat<-as.matrix(fit)
  stats<-apply(mat,2,function(xx)c('mean'=convertFunc(mean(xx)),'lower'=convertFunc(unname(quantile(xx,c(.025)))),'upper'=convertFunc(unname(quantile(xx,c(.975)))),'p'=findCredLim(xx)))
  return(stats)
}
findCredLim<-function(xx,target=0,tol=min(1/length(xx)/10,.001)){
  p<-optimize(function(quant)abs(quantile(xx,quant)-target),c(0,1),tol=tol)$minimum
  if(p>.5)p<-1-p
  p<-max(p*2,1/length(xx))
  return(p)
}

assignNames<-function(stats,virus){
  stats<-stats[,grep('foldChange|metaFoldChange\\[',colnames(stats))]
  foldIds<-grep('foldChange',colnames(stats))
  nums<-as.numeric(sub('.*\\[([0-9]+)\\]','\\1',colnames(stats)[foldIds]))
  colnames(stats)[foldIds]<-virus[nums]
  colnames(stats)[colnames(stats)=='metaFoldChange']<-'Overall'
  return(stats)
}

plotFit<-function(fit,virus,main='',cols=NULL,xlims=c(min(c(1,selects)),max(c(1,selects))),special=NULL,xlab='Fold change'){
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
  if(!is.null(special))selects<-selects[,order(colnames(selects) %in% special)]
  notInCols<-colnames(selects)[!colnames(selects) %in% names(cols)]
  cols<-c(cols,structure(rep('black',length(notInCols)),.Names=notInCols))
  plot(1,1,ylim=c(1,ncol(selects)),xlim=xlims,log='x',type='n',ylab='',xlab=xlab,main=main,yaxt='n',xaxt='n',mgp=c(2.5,.8,0),bty='l')
  axis(2,1:ncol(selects),colnames(selects),las=1)
  segments(selects['lower',],1:ncol(selects),selects['upper',],1:ncol(selects))
  points(selects['mean',],1:ncol(selects),pch=21,bg=cols[colnames(selects)])
  dnar::logAxis(1)
  if(xlims[2]<3&xlims[1]>.1)axis(1,.5)
  if(xlims[2]<10&xlims[1]>.5)axis(1,.7)
  lowerBound<-par('usr')[1]
  lows<-selects['lower',]<10^(lowerBound)
  if(any(lows)){
    for(ii in which(lows))polygon(c(10^(lowerBound),rep(10^(lowerBound+diff(par('usr')[1:2])*.04),2)),c(ii,ii+.1,ii-.1),col='black')
  }
  abline(v=1,lty=3)
  if(!is.null(special))abline(h=max(which(colnames(selects)%in% special))-.5)
}

compareAlleles<-function(mod,allele1,allele2,dat1,dat2=dat1,wtPar=NULL,logFunc=log,chains=50,...){
  message(allele1,' - ',allele2)
  dat=list(
    nEnv=nrow(dat1),
    nRep=3,
    prop1=logFunc(dat1[,grep(sprintf('^%s [0-9]+$',allele1),colnames(dat1))]/100),
    prop2=logFunc(dat2[,grep(sprintf('^%s [0-9]+$',allele2),colnames(dat2))]/100)
  )
  if(!is.null(wtPar))dat<-c(dat,list('wtMu'=wtPar[1],'wtSd'=wtPar[2]))
  fit <- sampling(mod, data = dat, iter=20000, chains=chains,thin=2,control=list(adapt_delta=.99,max_treedepth=15),...)
}


stanCodeSites<-'
  data {
    int<lower=0> nChimp;
    int<lower=0> nSites;
    int<lower=0> nAllele;
    int<lower=0,upper=1> siv[nChimp];
    //int<lower=0,upper=nSites> sites[nChimp];
    matrix<lower=0,upper=1>[nChimp,nAllele] alleles;
  }
  parameters {
    //vector[nSites] siteIntercepts;
    vector[nAllele] alleleBetas;
    //real<lower=0> shrinkage;
    real alpha;
    real<lower=0> tau;
    real<lower=0> lambda;
  }
  transformed parameters{
  }
  model {
    //siv~bernoulli_logit(alleles*alleleBetas+siteIntercepts[sites]);
    siv~bernoulli_logit(alpha+alleles*alleleBetas);
    //alleleBetas~double_exponential(0,shrinkage);
    //alleleBetas~normal(0,1);
    alleleBetas~normal(0,lambda*tau);
    lambda ~ cauchy(0,1);
    tau ~ cauchy(0,1);
    //siteIntercepts~normal(0,10);
    //shrinkage~inv_gamma(.001,.001);
  }
'
modSites <- stan_model(model_code = stanCodeSites)

analyzeSites<-function(siv,sites,alleles,chains=30,...){
  if(length(siv)!=length(sites)||length(siv)!=nrow(alleles))stop('Length differences in chimps and sites')
  siteNum<-as.numeric(as.factor(sites))
  dat<-list(
    nChimp=length(siv),
    nSites=max(siteNum),
    nAllele=ncol(alleles),
    siv=siv*1,
    sites=siteNum,
    alleles=alleles*1
  )
  fit <- sampling(modSites, data = dat, iter=100000, chains=chains,thin=4,control=list(adapt_delta=.95,max_treedepth=15),...)
}
