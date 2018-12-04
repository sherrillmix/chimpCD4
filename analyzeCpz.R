library('rstan')
envColors<-read.csv('envColor.csv',stringsAsFactors=FALSE)
envColors$col<-sprintf('#%s',envColors$color)
envCols<-structure(envColors$col,.Names=envColors$env)

back<-as.numeric(as.vector(read.csv('Background per allele for Scott.csv',stringsAsFactors=FALSE)[1,-1]))
#summary(glm(back~as.factor(rep(1:10,each=3))))
pdf('test.pdf');plot(1:length(back),back,xlab='Background sample',ylab='Background % infection');abline(v=seq(3,27,3)+.5);dev.off()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
cpz<-read.csv('SIVcpz Envs for Scott.csv',stringsAsFactors=FALSE)
colnames(cpz)[1]<-'env'
colnames(cpz)[2:ncol(cpz)]<-sprintf('%s %d',rep(colnames(cpz)[seq(2,ncol(cpz),3)],each=3),1:3)
cpz[cpz$env=='Gab2','env']<-'GAB2'
rownames(cpz)<-cpz$env
cpz<-cpz[dnar::orderIn(rownames(cpz),names(envCols)),]

cpzNorm<-cpz[,-1]/unlist(cpz[,2:4])


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
    real prop1Base[nEnv];
    vector[nEnv] metaRaw;
    real<lower=0> sigma;
    vector[nRep-1] repEffect;
    vector[nRep-1] prop1Noise[nEnv];
    real<lower=0> noiseSigma;
  }
  transformed parameters{
    vector[nEnv] foldChange;
    vector[nRep] repEffectPlusZero;
    vector[nRep] prop1NoisePlusZero[nEnv];
    foldChange=metaFoldChange+metaRaw*metaSd;
    repEffectPlusZero[1]=0;
    repEffectPlusZero[2:nRep]=repEffect[1:(nRep-1)];
    for(ii in 1:nEnv){
      prop1NoisePlusZero[ii][1]=0;
      prop1NoisePlusZero[ii][2:nRep]=prop1Noise[ii];
    }
  }
  model {
    for(ii in 1:nEnv){
      prop1Noise[ii,]~normal(0,noiseSigma);
      prop1[ii,]~normal(prop1Base[ii]+repEffectPlusZero+prop1NoisePlusZero[ii],sigma);
      prop2[ii,]~normal(prop1Base[ii]+foldChange[ii]+repEffectPlusZero+prop1NoisePlusZero[ii],sigma);
    }
    //back~normal(-1.8,1);
    metaRaw~normal(0,1);
  }
'
mod <- stan_model(model_code = stanCode)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
uniqAlleles<-unique(sub(' .*$','',colnames(cpz)[-1]))
names(uniqAlleles)<-uniqAlleles
pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('Human','QQKVP'),c('QQKVP','QQKVT'),c('QQNVP','RQNVP'),c('QQNVP','QRNVP'),c('QQNIP','QRNIP'),c('QQKVT','QRKVT'),c('QQNVP','QQKVP'),c('QQNVT','QQKVT'),c('QQNVP','QQNIP'),c('QRNVP','QRNIP'))
compareAlleles<-function(allele1,allele2){
  message(allele1,' - ',allele2)
  dat=list(
    nEnv=nrow(cpz),
    nRep=3,
    prop1=log(cpz[,grep(allele1,colnames(cpz))]/100),
    prop2=log(cpz[,grep(allele2,colnames(cpz))]/100)
  )
  fit <- sampling(mod, data = dat, iter=10000, chains=10,thin=20,control=list(adapt_delta=.9))
}

plotRaw<-function(raw1,raw2,lab1=sub(' .*$','',colnames(raw1)[1]),lab2=sub(' .*$','',colnames(raw2)[1]),ylim=range(cbind(raw1,raw2)),cols=NULL){
  par(mfrow=c(2,ceiling(nrow(raw1)/2))) 
  notInCols<-rownames(raw1)[!rownames(raw1) %in% names(cols)]
  cols<-c(cols,structure(rep('black',length(notInCols)),.Names=notInCols))
  for(ii in 1:nrow(raw1)){
    plot(1,1,type='n',xlim=c(.5,2.5),ylim=ylim,log='y',yaxt='n',xlab='',ylab='Percent infected',main=rownames(raw1)[ii],xaxt='n')
    axis(1,1:2,c(lab1,lab2))
    points(rep(1,ncol(raw1)),raw1[ii,],col=cols[ii])
    points(rep(2,ncol(raw2)),raw2[ii,],col=cols[ii])
    segments(rep(1,ncol(raw1)),unlist(raw1[ii,]),rep(2,ncol(raw2)),unlist(raw2[ii,]))
    text(rep(1.5,ncol(raw1)),exp((log(raw1[ii,])+log(raw2[ii,]))/2),1:ncol(raw1))
    dnar::logAxis(las=1)
  }
}

plotFit<-function(fit,virus,main='',cols=NULL){
  mat<-as.matrix(fit)
  stats<-apply(mat,2,function(xx)c('mean'=mean(xx),'lower'=unname(quantile(xx,c(.025))),'upper'=unname(quantile(xx,c(.975)))))
  selects<-exp(stats[,grep('foldChange|metaFoldChange|repEffect\\[',colnames(stats))])
  foldIds<-grep('foldChange',colnames(selects))
  nums<-as.numeric(sub('.*\\[([0-9]+)\\]','\\1',colnames(selects)[foldIds]))
  colnames(selects)[foldIds]<-virus[nums]
  colnames(selects)[colnames(selects)=='metaFoldChange']<-'Overall'
  colnames(selects)[colnames(selects)=='repEffect[1]']<-'Replicate 2'
  colnames(selects)[colnames(selects)=='repEffect[2]']<-'Replicate 3'
  selects<-selects[,ncol(selects):1]
  notInCols<-colnames(selects)[!colnames(selects) %in% names(cols)]
  cols<-c(cols,structure(rep('black',length(notInCols)),.Names=notInCols))
  par(mar=c(3.5,6,2,.4))
  xlims<-c(min(c(1,selects)),max(c(1,selects)))
  plot(1,1,ylim=c(1,ncol(selects)),xlim=xlims,log='x',type='n',ylab='',xlab='Fold change',main=main,yaxt='n',xaxt='n',mgp=c(2.5,.8,0))
  axis(2,1:ncol(selects),colnames(selects),las=1)
  segments(selects['lower',],1:ncol(selects),selects['upper',],1:ncol(selects))
  points(selects['mean',],1:ncol(selects),pch=21,bg=cols[colnames(selects)])
  dnar::logAxis(1)
  if(xlims[2]<3&xlims[1]>.1)axis(1,.5)
  if(xlims[2]<10&xlims[1]>.5)axis(1,.7)
  abline(v=1,lty=3)
}

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair)compareAlleles(pair[1],pair[2]))
  names(selectFits)<-sapply(pairs,paste,collapse='-')
}

dummy<-lapply(names(selectFits),function(xx){message(xx);print(selectFits[[xx]],pars=c('metaFoldChange','foldChange','repEffect'))})
pdf('cpzFits.pdf',width=4,height=4);lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(cpz),main=xx,cols=envCols));dev.off()
pdf('cpzRaw.pdf',width=8,height=6);lapply(pairs,function(xx)plotRaw(cpz[,grep(xx[1],colnames(cpz))],cpz[,grep(xx[2],colnames(cpz))],cols=envCols));dev.off()

#allFits<-lapply(uniqAlleles,function(allele1){
  #lapply(uniqAlleles,compareAlleles,allele1)
#})


#print(fit,pars=c('metaFoldChange','foldChange','repEffect'))
#print(plot(fit))
#library('bayesplot')
#mcmc_areas(as.matrix(fit),regex_pars=c('foldChange.*','repEffect.*'),prob=.95)
