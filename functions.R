library('rstan')
if(!dir.exists('out'))dir.create('out')
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

envColors<-read.csv('data/envColor.csv',stringsAsFactors=FALSE)
envColors$col<-sprintf('#%s',envColors$color)
envCols<-structure(envColors$col,.Names=envColors$env)


stanCodeCpz<-'
  data {
    int<lower=0> nEnv;
    int<lower=0> nRep;
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
modCpz <- stan_model(model_code = stanCodeCpz)

stanCodeMonkey<-'
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
modMonkey <- stan_model(model_code = stanCodeMonkey)

stanCodeMB<-'
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
    }
    for(ii in 1:nEnv){
      back1[ii,]~normal(-6.43,.97);
      back2[ii,]~normal(-6.43,.97);
    }
    metaRaw~normal(0,1);
    wtFoldChange~normal(wtMu,wtSd);
  }
'
modMB <- stan_model(model_code = stanCodeMB)

logAxis<-function(side=2,exponent=TRUE,addExtra=!exponent,minorTcl=-.2,axisMin=-Inf,axisMax=Inf,offset=0,col.ticks='black',...){
  if(side %in% c(2,4)) parX<-sort(graphics::par('usr')[3:4])
  else parX<-sort(graphics::par('usr')[1:2])
  minX<-max(10^parX[1],axisMin)
  maxX<-min(10^parX[2],axisMax)
  if(log10(maxX)-log10(minX)>400)stop(simpleError('Huge range in logged axis'))
  allTicks<-unlist(lapply(floor(log10(minX)):ceiling(log10(maxX)),function(x)1:9*10^x))
  allTicks<-allTicks[allTicks<=maxX & allTicks>=minX]
  graphics::axis(side,allTicks+offset,rep('',length(allTicks)),tcl=minorTcl,col.ticks=col.ticks)
  if(ceiling(log10(minX))<=floor(log10(maxX)))prettyY<-seq(ceiling(log10(minX)),floor(log10(maxX)),1)
  else prettyY<-c()
  graphics::axis(side,10^prettyY+offset,rep('',length(prettyY)),tcl=minorTcl*2,col.ticks=col.ticks)
  if(length(prettyY)>7)prettyY<-pretty(prettyY)
  if(length(prettyY)==0)prettyY<-c(ceiling(log10(minX)),floor(log10(maxX)))
  if(addExtra){
    origPretty<-prettyY
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(5),origPretty-log10(10/5)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(2),origPretty-log10(10/2)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(3),origPretty-log10(10/3)))
    if(sum(prettyY>=log10(minX)&prettyY<=log10(maxX))<4)prettyY<-unique(c(prettyY,origPretty+log10(7),origPretty-log10(10/7)))
  }
  if(exponent){
    if(any(prettyY%%1!=0))labs<-sapply(prettyY,function(x)as.expression(bquote(.(10^(x%%1))%*%10^.(floor(x)))))
    else labs<-ifelse(prettyY==0,1,sapply(prettyY,function(x)as.expression(bquote(10^.(floor(x))))))
  }
  else labs<-10^prettyY
  graphics::axis(side,10^prettyY+offset,labs,col.ticks=col.ticks,...)
  return(invisible(list('minor'=allTicks,'major'=10^prettyY)))
}

plotRaw<-function(raw1,raw2,lab1=sub(' .*$','',colnames(raw1)[1]),lab2=sub(' .*$','',colnames(raw2)[1]),ylim=range(cbind(raw1,raw2,10)),cols=NULL){
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
    logAxis(las=1)
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
  logAxis(1)
  if(xlims[2]<3&xlims[1]>.1)axis(1,.5)
  if(xlims[2]<10&xlims[1]>.5)axis(1,.7)
  lowerBound<-par('usr')[1]
  lows<-selects['lower',]<10^(lowerBound)
  if(any(lows)){
    for(ii in which(lows))polygon(c(10^(lowerBound),rep(10^(lowerBound+diff(par('usr')[1:2])*.04),2)),c(ii,ii+.1,ii-.1),col='black')
  }
  abline(v=1,lty=3)
  if(!is.null(special))abline(h=max(which(colnames(selects)%in% special))-.5,lty=2)
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

plotGlm<-function(mod,xlim=range(c(mod$upper,mod$lower)),main='',xlab='Fold change in odds of SIVcpz infection',cex.axis=.9,...){
  mod$upper<-mod$Estimate+mod[,'Std. Error']*2
  mod$lower<-mod$Estimate-mod[,'Std. Error']*2
  ylims<-c(.5,nrow(mod))
  plot(1,1,type='n',xlim=exp(xlim),ylim=ylims,ylab='',xlab=xlab,yaxt='n',log='x',xaxt='n',mgp=c(2.25,1,0),bty='l',...)
  par(lheight=.7)
  mtext(sub(' ','\n',main),3,xpd=NA,at=ifelse(nchar(main)>20,1,1),line=ifelse(grepl(' ',main),0,.7),cex=cex.axis/1.2)
  points(exp(mod$Estimate),nrow(mod):1,pch=21,bg='black',cex=.8)
  segments(exp(mod$lower),nrow(mod):1,exp(mod$upper),nrow(mod):1)
  axis(2,nrow(mod):1,rownames(mod),las=1,xpd=NA,tcl=-.2,mgp=c(3,.5,0),cex.axis=cex.axis)
  logAxis(1,mgp=c(3,.4,0),cex.axis=cex.axis)
  abline(v=1,lty=2)
}

