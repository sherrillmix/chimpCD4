set.seed(12345)
source('functions.R')

alleles<-dnar::fillDown(unlist(read.csv('dat/Mus CD4 data for Scott.csv',nrows=1,header=FALSE,stringsAsFactors=FALSE,check.names=FALSE))[-1])
dates<-unlist(read.csv('dat/Mus CD4 data for Scott.csv',nrows=2,header=FALSE,stringsAsFactors=FALSE)[2,])[-1]
dat<-read.csv('dat/Mus CD4 data for Scott.csv',stringsAsFactors=FALSE,row.names=1,skip=2,header=FALSE)
monk<-data.frame('rep'=unlist(dat),'virus'=rep(rownames(dat),ncol(dat)),'date'=rep(dates,each=nrow(dat)),'allele'=rep(alleles,each=nrow(dat)),stringsAsFactors=FALSE)
monk<-monk[monk$virus!='SIVden CD1',]
monk$siv<-ifelse(!grepl('SIV',monk$virus),'SIVcpz',ifelse(grepl('SIVgor',monk$virus),'SIVgor','Monkey SIV'))
monk$siv[grepl('Mock',monk$virus)]<-'Mock'
monk$species<-ifelse(monk$siv=='SIVcpz','cpz',ifelse(monk$siv=='SIVgor','gor',sub('SIV(...).*','\\1',monk$virus)))
monk$repNo0<-ifelse(monk$rep==0&!is.na(monk$rep),min(monk$rep[monk$rep>0],na.rm=TRUE),monk$rep)

bono<-read.csv('dat/bonobo CD4 data for Scott.csv')
alleles<-dnar::fillDown(unlist(read.csv('dat/gorilla CD4 data for Scott 10162020.csv',nrows=1,header=FALSE,stringsAsFactors=FALSE)))[-1]
dates<-unlist(read.csv('dat/gorilla CD4 data for Scott 10162020.csv',nrows=2,header=FALSE,stringsAsFactors=FALSE)[2,])[-1]
dat<-read.csv('dat/gorilla CD4 data for Scott 10162020.csv',row.names=1,skip=2,header=FALSE,stringsAsFactors=FALSE)
gor<-data.frame('rep'=unlist(dat),'virus'=rep(rownames(dat),ncol(dat)),'date'=rep(dates,each=nrow(dat)),'allele'=rep(alleles,each=nrow(dat)),stringsAsFactors=FALSE)
gor$siv<-ifelse(!grepl('SIV',gor$virus),'SIVcpz',ifelse(grepl('SIVgor',gor$virus),'SIVgor','Monkey SIV'))
gor$siv[grepl('Mock',gor$virus)]<-'Mock'
gor$species<-ifelse(gor$siv=='SIVcpz','cpz',ifelse(gor$siv=='SIVgor','gor',sub('SIV(...).*','\\1',gor$virus)))
hum<-gor[gor$allele=='Human',]
rownames(hum)<-paste(hum$virus,hum$date)
gor$hum<-hum[paste(gor$virus,gor$date),'rep']
gor$norm<-gor$rep/gor$hum

means<-tapply(gor$norm,gor[colnames(gor)!='Human',c('virus','allele')],mean,na.rm=TRUE)
apply(means[,!grepl('Bonobo|Human',colnames(means))],2,function(xx){(wilcox.test(xx[unique(gor[gor$siv=='SIVcpz','virus'])],xx[unique(gor[gor$siv=='SIVgor','virus'])]))$p.value})
apply(means[,!grepl('Bonobo|Human',colnames(means))],2,function(xx){(t.test(xx[unique(gor[gor$siv=='Monkey SIV','virus'])],xx[unique(gor[gor$siv=='SIVgor','virus'])]))$p.value})
apply(means[,!grepl('Bonobo|Human|N15T',colnames(means))],2,function(xx){(wilcox.test(xx[unique(gor[gor$siv=='Monkey SIV','virus'])],xx[unique(gor[gor$siv=='SIVcpz','virus'])]))$p.value})
apply(means[,!grepl('Bonobo|Human|N15T',colnames(means))],2,function(xx){(var.test(xx[unique(gor[gor$siv=='Monkey SIV','virus'])],xx[unique(gor[gor$siv=='SIVcpz','virus'])]))$p.value})
sds<-apply(means[!grepl('Mock',rownames(means)),!grepl('Bonobo|Human|N15T',colnames(means))],1,function(xx)sd((xx)))
#sdsLog<-apply(means[!grepl('Mock',rownames(means)),!grepl('Bonobo|Human|N15T',colnames(means))],1,function(xx)sd(log(xx)))
#geoSd<-function(xx)exp(sqrt(mean(log(xx/exp(mean(log(xx))))^2)))
geoSd<-function(xx)exp(sd(log(xx)))
sdsLog<-apply(means[!grepl('Mock',rownames(means)),!grepl('Bonobo|Human|N15T',colnames(means))],1,function(xx)geoSd(xx))
t.test(sds[unique(gor[gor$siv=='Monkey SIV','virus'])],sds[unique(gor[gor$siv=='SIVcpz','virus'])])
t.test(sds[unique(gor[gor$siv=='SIVgor','virus'])],sds[unique(gor[gor$siv=='SIVcpz','virus'])])
t.test(sds[unique(gor[gor$siv=='Monkey SIV','virus'])],sds[unique(gor[gor$siv=='SIVgor','virus'])])
t.test(sdsLog[unique(gor[gor$siv=='Monkey SIV','virus'])],sdsLog[unique(gor[gor$siv=='SIVcpz','virus'])])
t.test(sdsLog[unique(gor[gor$siv=='SIVgor','virus'])],sdsLog[unique(gor[gor$siv=='SIVcpz','virus'])])
t.test(sdsLog[unique(gor[gor$siv=='Monkey SIV','virus'])],sdsLog[unique(gor[gor$siv=='SIVgor','virus'])])
wilcox.test(sdsLog[unique(gor[gor$siv=='Monkey SIV','virus'])],sdsLog[unique(gor[gor$siv=='SIVcpz','virus'])])
wilcox.test(sdsLog[unique(gor[gor$siv=='Monkey SIV','virus'])],sdsLog[unique(gor[gor$siv=='SIVgor','virus'])])
#wilcox.test(sdsLog[unique(gor[gor$siv=='SIVgorMonkey SIV','virus'])],sdsLog[unique(gor[gor$siv=='SIVcpz','virus'])])
pos<-structure(1:3,.Names=c('SIVgor','SIVcpz','Monkey SIV'))
xPos<-pos[sapply(names(sds),function(xx)gor[gor$virus==xx,'siv'][1])]
pdf('sd.pdf')
plot(xPos,sds,xlim=c(.5,3.5),xaxt='n',ylab='Standard deviation',xlab='',ylim=c(0,max(sds)),las=1)
text(xPos+.02,sds,names(sds),cex=.5,adj=0)
axis(1,xPos,names(xPos))
plot(xPos,sdsLog,xlim=c(.5,3.5),xaxt='n',ylab='Geometric standard deviation',xlab='',ylim=c(0,max(sdsLog)),las=1)
text(xPos+.02,sdsLog,names(sds),cex=.5,adj=0)
axis(1,xPos,names(xPos))
dev.off()

#https://mokole.com/palette.html
cols<-c('#2f4f4f', '#2e8b57', '#8b0000', '#808000', '#000080', '#ff0000', '#ff8c00', '#ffd700', '#c71585', '#00ff7f', '#0000ff', '#f08080', '#adff2f', '#ff00ff', '#1e90ff', '#dda0dd', '#87ceeb', '#7fffd4', '#ffe4c4','black')
virCols<-structure(cols,.Names=unique(gor$virus))
allelePairs<-list(c('AHSM','ARSM','ANSR'),c('PRSM','AHPM'))
allelePos<-structure(1:length(unlist(allelePairs)),.Names=unlist(allelePairs))
pch=structure(c(20+1:length(unique(gor$siv))),.Names=unique(gor$siv))
pdf('testGor.pdf',height=3,width=8)
par(mfrow=c(1,3),mar=c(8,4,1.2,.8))
for(ii in unique(gor$siv[!grepl('Mock',gor$siv)])){
  plot(1,1,type='n',bty='l',main=ii,xlim=c(.5,max(allelePos)+.5),xaxs='i',xaxt='n',ylim=c(0,100),xlab='',ylab='Infectivity (% human allele)',xpd=NA,las=1,mgp=c(2.2,.7,0))
  dnar::slantAxis(1,allelePos,names(allelePos),location=1,srt=-30)
  thisMeans<-means[rownames(means) %in% gor[gor$siv==ii,'virus'],,drop=FALSE]
  for(jj in allelePairs){
    sapply(rownames(thisMeans),function(xx)lines(allelePos[jj],thisMeans[xx,jj]*100,col=virCols[xx]))
  }
  thisCols<-virCols[names(virCols) %in% rownames(thisMeans)]
  points(rep(allelePos,each=nrow(thisMeans)),unlist(thisMeans[,names(allelePos)])*100,pch=pch[ii],bg=thisCols[rep(rownames(thisMeans),ncol(thisMeans))],cex=1.2)
  legend(mean(par('usr')[1:2]),dnar::convertLineToUser(3.3),names(thisCols),pch=21,xpd=NA,pt.bg=thisCols,pt.cex=1.2,ncol=2,yjust=1,xjust=.5,bty='n')
}
dev.off()
pdf('test.pdf')
datePos<-structure(1:length(unique(gor$date)),.Names=unique(gor$date))
plot(datePos[gor[gor$allele=='Human','date']],gor[gor$allele=='Human','rep'])
tapply(gor[gor$allele=='Human','rep'],gor[gor$allele=='Human',c('virus')],function(xx)lines(datePos,xx)) #lazy
dev.off()

pdf('heat.pdf');heatmap(means[!grepl('Mock',rownames(means)),],scale='none');dev.off()
#rep~virus*allele+species+virus+allele+date


stanCodeGor<-'
  data {
    int<lower=0> nRep;
    vector[nRep] prop;
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
modGor <- stan_model(model_code = stanCodeGor)


stanCodeGor<-'
  data {
    int<lower=0> nObs;
    int<lower=0> nHuman;
    int<lower=0> nVirus;
    int<lower=0> nAllele;
    vector[nHuman]<lower=0> human;
    //rep data for each virus-day
    int<lower=0> virusId[nObs];
    int<lower=0> humanDateId[nObs];
    vector[nObs]<lower=0> reps;
    //human data for each virus-day 
    int<lower=0> humanReps[nHuman];
    int<lower=0> dateId[nHuman];
    int<lower=0> humanVirusId[nHuman];
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

virLookup<-list(
  'SIVgor'=c('SIVgorCP2135'='CP2135','SIVgorCP2139'='CP2139','SIVgorBQID2'='BQID2','SIVgorBPID1'='BPID1'),
  'SIVcpz'=c('MT145'='MT145','MB897'='MB897','EK505'='EK505','LB715'='LB715','GAB2'='GAB2','TAN2'='TAN2','BF1167'='BF1167'),
  'Other'=c('SIVmus1-54'='SIVmus1-54','SIVmusGAB11'='SIVmusGAB11','SIVlst7'='SIVlstUS7','SIVagmTAN1'='SIVagmTAN1','SIVwrcGM05'='SIVwrcGM05','SIVasc'='SIVascRT11','SIVsmmFTq'='SIVsmmFTq','SIVsmmSL92b'='SIVsmmSL92b')
)
alleleLookup=c('Human','AHSM','ARSM','ANSR','PRSM','AHPM')

pdf('heat.pdf')
tmp<-means[!grepl('Mock',rownames(means)),!grepl('Human',colnames(means))]
method<-'euclidean'
virClust<-hclust(dist((tmp),method=method))
alleleClust<-hclust(dist(t((tmp)),method=method))
virClust2<-hclust(dist((tmp[,alleleLookup[alleleLookup!='Human']]),method=method))
tmp<-tmp[virClust$labels[virClust$order],alleleClust$labels[alleleClust$order]]
#heatmap(tmp,scale='none');
plotFunc<-function(tmp){
  breaks<-seq(0,1,length.out=101)
  #cols<-hcl.colors(100, "YlOrRd", rev = TRUE)
  cols<-colorRampPalette(c('white','#f7c252','#eb5500','#7d0025'))(100)
  image(1:ncol(tmp),1:nrow(tmp),t(tmp),xaxt='n',yaxt='n',ylab='',xlab='',breaks=breaks,col=cols)
  dnar::slantAxis(1,1:ncol(tmp),colnames(tmp))
  axis(2,1:nrow(tmp),rownames(tmp),las=1)
  par(lheight=.7)
  dnar::insetScale(breaks,cols,main='Proportion of\nhuman replication')
  abline(h=2:nrow(means)-.5,v=2:nrow(means)-.5,col='#FFFFFF33',lwd=2)
  box()
}

par(mar=c(5,7,.1,.3))
plotFunc(tmp)
virOrder<-rev(unlist(lapply(virLookup,names)))
tmp2<-tmp[virOrder,alleleLookup[alleleLookup!='Human']]
rownames(tmp2)<-virOrder
plotFunc(tmp2)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
tmp<-tmp[virClust2$labels[virClust$order],alleleClust$labels[alleleClust$order]]
tmp3<-tmp[order(sapply(rownames(tmp),function(xx)-which(sapply(virLookup,function(zz)xx%in%names(zz))))),colnames(tmp)[colnames(tmp) %in% colnames(tmp2)]]
#tmp3<-tmp[order(sapply(rownames(tmp),function(xx)-which(sapply(virLookup,function(zz)xx%in%names(zz)))),-apply(tmp,1,sum)),colnames(tmp)[colnames(tmp) %in% colnames(tmp2)]]
plotFunc(tmp3)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
tmp4<-tmp[order(sapply(rownames(tmp),function(xx)-which(sapply(virLookup,function(zz)xx%in%names(zz))))),c(colnames(tmp)[colnames(tmp) %in% colnames(tmp2)],'AHSM N15T')]
par(mar=c(5,7,.1,2))
plotFunc(tmp4)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
abline(v=5.5,lwd=2)
dev.off()


virLookup<-list(
  'Other'=c('SIVmusGAB11'='SIVmusGAB11','SIVmus1-54'='SIVmus1-54','SIVgsn CN71'='SIVgsnCN71','SIVmon'='SIVmonCM99','SIVdeb'='SIVdeb04CMP','SIVlst7'='SIVlstUS7','SIVagmTAN1'='SIVagmTAN1','SIVwrcGM05'='SIVwrcGM05','SIVasc'='SIVasc','SIVsmmFTq'='SIVsmmFTq','SIVsmmSL92b'='SIVsmmSL92b'),
  'SIVcpz'=c('MT145'='SIVcpzMT145','MB897'='SIVcpzMB897','BF1167'='SIVcpzBF1167','TAN2'='SIVcpzTAN2'),
  'SIVgor'=c('SIVgorCP2135'='SIVgorCP2135','SIVgorBPID1'='SIVgorBPID1','SIVgorBQID2'='SIVgorBQID2')
)
alleleLookup<-unique(monk$allele)
means<-tapply(monk$rep,monk[,c('virus','allele')],mean,na.rm=TRUE)

pdf('Fig._3_heat.pdf')
tmp<-means[!grepl('Mock',rownames(means)),]
method<-'euclidean'
virClust<-hclust(dist(log(tmp),method=method),method='average')
#weights<-ifelse((virClust$labels[virClust$order]) %in% c('SIVmus1-54','SIVmusGAB11'),500,1)
weights<-ifelse(virClust$labels=='SIVmusGAB11',1000,ifelse(virClust$labels=='SIVmus1-54',500,1))
library(vegan)
virClust<-reorder(virClust,weights)
alleleClust<-hclust(dist(t(log(tmp)),method=method),method='average')
weights<-ifelse(alleleClust$labels %in% c('Human'),-50,1)
alleleClust<-reorder(alleleClust,weights)
tmp<-tmp[virClust$labels[virClust$order],alleleClust$labels[alleleClust$order]]
#heatmap(tmp,scale='none');
plotFunc<-function(tmp){
  breaks<-seq(0,max(means),length.out=101)
  #cols<-hcl.colors(100, "YlOrRd", rev = TRUE)
  cols<-colorRampPalette(c('white','#f7c252','#eb5500','#7d0025'))(100)
  image(1:ncol(tmp),1:nrow(tmp),t(tmp),xaxt='n',yaxt='n',ylab='',xlab='',breaks=breaks,col=cols)
  dnar::slantAxis(1,1:ncol(tmp),colnames(tmp))
  axis(2,1:nrow(tmp),rownames(tmp),las=1)
  par(lheight=.7)
  dnar::insetScale(breaks,cols,main='Percentage of\ncells infected')
  abline(h=2:nrow(means)-.5,v=2:nrow(means)-.5,col='#FFFFFF33',lwd=2)
  box()
}
par(mar=c(5,7,.1,.6))
plotFunc(tmp)
virOrder<-rev(unlist(lapply(virLookup,names)))
tmp2<-tmp[virOrder,alleleLookup]
rownames(tmp2)<-do.call(c,unname(virLookup))[rownames(tmp2)]
plotFunc(tmp2)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
virOrder<-virClust$labels[virClust$order]
virOrder<-virOrder[virOrder!='SIVden CD1']
tmp<-tmp[virOrder,alleleClust$labels[alleleClust$order]]
tmp3<-tmp[order(sapply(rownames(tmp),function(xx)-which(sapply(virLookup,function(zz)xx%in%names(zz))))),colnames(tmp)[colnames(tmp) %in% colnames(tmp2)]]
#tmp3<-tmp[order(sapply(rownames(tmp),function(xx)-which(sapply(virLookup,function(zz)xx%in%names(zz)))),-apply(tmp,1,sum)),colnames(tmp)[colnames(tmp) %in% colnames(tmp2)]]
rownames(tmp3)<-do.call(c,unname(virLookup))[rownames(tmp3)]
plotFunc(tmp3)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
virOrder<-rev(unlist(lapply(virLookup,names)))
tmp<-tmp[virOrder,alleleClust$labels[alleleClust$order]]
tmp4<-tmp[virOrder,colnames(tmp)[colnames(tmp) %in% colnames(tmp2)]]
rownames(tmp4)<-do.call(c,unname(virLookup))[rownames(tmp4)]
plotFunc(tmp4)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
dev.off()


sort(sapply(structure(unique(gor$virus),.Names=unique(gor$virus)),function(ii)t.test(log(gor[gor$allele=='AHSM'&gor$virus==ii,'rep']),log(gor[gor$allele=='AHSM N15T'&gor$virus==ii,'rep']))$p.value))
sort(sapply(structure(unique(gor$virus),.Names=unique(gor$virus)),function(ii)t.test(gor[gor$allele=='Bonobo I83'&gor$virus==ii,'rep'],gor[gor$allele=='Bonobo T83'&gor$virus==ii,'rep'],paired=TRUE)$p.value))
t.test(log(gor[gor$allele=='Bonobo I83'&!grepl('Mock',gor$virus),'rep']),log(gor[gor$allele=='Bonobo T83'&!grepl('Mock',gor$virus),'rep']),paired=TRUE)
wilcox.test(gor[gor$allele=='Bonobo I83'&!grepl('Mock',gor$virus),'rep'],gor[gor$allele=='Bonobo T83'&!grepl('Mock',gor$virus),'rep'],paired=TRUE)
wilcox.test(
  dnar::withAs(xx=gor[gor$allele=='Bonobo I83'&!grepl('Mock',gor$virus),],tapply(xx$rep,xx$virus,mean,na.rm=TRUE)),
  dnar::withAs(xx=gor[gor$allele=='Bonobo T83'&!grepl('Mock',gor$virus),],tapply(xx$rep,xx$virus,mean,na.rm=TRUE)),
paired=TRUE)
wilcox.test(dnar::withAs(xx=gor[gor$allele=='Bonobo T83'&!grepl('Mock',gor$virus),],tapply(xx$rep,xx$virus,mean,na.rm=TRUE)))

virs<-unique(gor$virus[!grepl('Mock',gor$virus)])
virs<-structure(virs,.Names=virs)
bonoPs<-do.call(rbind,lapply(virs,function(xx)c(
    'HumvsT83'=t.test(gor[gor$allele=='Human'&gor$virus==xx,'rep'],gor[gor$allele=='Bonobo T83'&gor$virus==xx,'rep'])$p.value,
    'HumvsI83'=t.test(gor[gor$allele=='Human'&gor$virus==xx,'rep'],gor[gor$allele=='Bonobo I83'&gor$virus==xx,'rep'])$p.value,
    'I83vsT83'=t.test(gor[gor$allele=='Bonobo I83'&gor$virus==xx,'rep'],gor[gor$allele=='Bonobo T83'&gor$virus==xx,'rep'])$p.value
  )
))
alls<-unique(gor$allele[!gor$allele %in% c('Bonobo T83','Bonobo I83','AHSM N15T')])
alls<-structure(alls,.Names=alls)
gorPs<-lapply(virs,function(vv)outer(alls,alls,function(xx,yy)mapply(function(xxx,yyy)t.test(gor[gor$allele==xxx&gor$virus==vv,'rep'],gor[gor$allele==yyy&gor$virus==vv,'rep'])$p.value,xx,yy)))
gorPStack<-do.call(rbind,lapply(names(gorPs),function(xx)data.frame('virus'=xx,'allele1'=matrix(rownames(gorPs[[xx]]),ncol=ncol(gorPs[[xx]]),nrow=nrow(gorPs[[xx]]))[upper.tri(gorPs[[xx]])],'allele2'=matrix(colnames(gorPs[[xx]]),ncol=ncol(gorPs[[xx]]),nrow=nrow(gorPs[[xx]]),byrow=TRUE)[upper.tri(gorPs[[xx]])],'p'=gorPs[[xx]][upper.tri(gorPs[[xx]])])))
write.csv(gorPStack,'allGorP.csv',row.names=FALSE)


t.test(
  dnar::withAs(xx=gor[gor$allele=='AHSM'&!grepl('Mock',gor$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE)),
  dnar::withAs(xx=gor[gor$allele=='AHSM N15T'&!grepl('Mock',gor$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE)),
paired=TRUE)
exp(-mean(dnar::withAs(xx=gor[gor$allele=='AHSM'&!grepl('Mock',gor$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE))-dnar::withAs(xx=gor[gor$allele=='AHSM N15T'&!grepl('Mock',gor$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE))))

combo<-cbind('i83'=gor[gor$allele=='Bonobo I83'&!grepl('Mock',gor$virus),c('virus','rep','date')],'t83'=gor[gor$allele=='Bonobo T83'&!grepl('Mock',gor$virus),c('virus','rep','date','siv')])
if(any(combo$i83.date!=combo$t83.date)||any(combo$i83.virus!=combo$t83.virus))stop('Mismatch')
combo$diff<-combo$t83.rep-combo$i83.rep



wilcox.test(
  dnar::withAs(xx=monk[monk$allele=='Allele 1'&!grepl('Mock',monk$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE)),
  dnar::withAs(xx=monk[monk$allele=='Allele 3'&!grepl('Mock',monk$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE)),
paired=TRUE)
wilcox.test(
  dnar::withAs(xx=monk[monk$allele=='Allele 2'&!grepl('Mock',monk$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE)),
  dnar::withAs(xx=monk[monk$allele=='Allele 4'&!grepl('Mock',monk$virus),],tapply(log(xx$rep),xx$virus,mean,na.rm=TRUE)),
paired=TRUE)


stanCode<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nAllele;
    int<lower=0> nRep;
    int<lower=0> nGroup;
    int<lower=0> nExperiment;
    real rep[nRep];
    int<lower=0> allele[nRep];
    int<lower=0> virus[nRep];
    int<lower=0> virusGroup[nVirus];
    int<lower=0> experiment[nRep];
    real backMean;
    real<lower=0> backSd;
    //int<lower=0> nMock;
    //real<lower=0> mock[nMock];
  }
  parameters {
    matrix[nAllele,nVirus] alleleVirus;
    matrix[nVirus,nExperiment-1] virusExperiment;
    matrix[nAllele-1,nGroup] alleleGroup;
    real<lower=0> sdGroup[nGroup];
    real<lower=0> repSd;
    vector<lower=0>[nRep] background;
    //real<lower=0> humMean;
    //real<lower=0> humSd;
    real<lower=0> alleleSd;
    vector[nAllele] alleleMeans;
    //real<lower=0> nu;
    //real<lower=0> backMean;
    //real<lower=0> backSd;
  }
  transformed parameters{
    vector[nRep] beta;
    for(ii in 1:nRep){
      beta[ii]=alleleVirus[1,virus[ii]];
      if(experiment[ii]!=1)beta[ii]=beta[ii]+virusExperiment[virus[ii],experiment[ii]-1];
      if(allele[ii]!=1)beta[ii]=beta[ii]+alleleVirus[allele[ii],virus[ii]];
      //beta[ii]=log_sum_exp(beta[ii],log(background)[ii]);
    }
    //beta=log(exp(beta)+background);
  }
  model {
    rep~normal(log(exp(beta)+background),repSd);
    for(ii in 1:nVirus)alleleVirus[2:nAllele,ii]~normal(alleleGroup[,virusGroup[ii]],sdGroup[virusGroup[ii]]);
    //for(ii in 1:nVirus)alleleVirus[2:nAllele,ii]~student_t(10,alleleGroup[,virusGroup[ii]],sdGroup[virusGroup[ii]]);
    for(ii in 1:nGroup)alleleGroup[,ii]~normal(alleleMeans,alleleSd);
    for(ii in 2:nExperiment)virusExperiment[,ii-1]~normal(0,10);
    repSd~gamma(1,.1);
    sdGroup~gamma(1,.1);
    background~normal(backMean,backSd);
    alleleVirus[1,]~normal(-2,10);
    alleleMeans~normal(0,10);
    alleleSd~gamma(1,.1);
    //nu~gamma(2,.1);
  }
'
modAllele <- stan_model(model_code = stanCode)

#mocks<-c(monk[monk$virus=='Mock'&!is.na(monk$rep),'rep'])
mocks<-c(gor[grepl('Mock',gor$virus)&!is.na(gor$rep),'rep'],monk[monk$virus=='Mock'&!is.na(monk$rep),'rep'])
fit<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$siv,xx$rep,xx$date,mocks,modAllele,chains=50))

save(fit,file='tmp.Rdat')
print(fit$fit,'alleleMeans')
print(fit$fit,c('beta','background'),include=FALSE)

mat<-as.matrix(fit$fit)
pdf('test.pdf');plot(monk[monk$virus!='Mock'&!is.na(monk$rep),'rep'],exp(apply(mat[,grep('beta',colnames(mat))],2,mean)),log='xy');dev.off()
diffs<-apply(mat[,!grepl('\\[1,',colnames(mat))&grepl('alleleVirus',colnames(mat))],2,quantile,c(.025,.975))
pdf('test.pdf')
plot(1,1,xlim=c(1,ncol(diffs)),ylim=range(diffs),type='n')
segments(1:ncol(diffs),diffs[1,],1:ncol(diffs),diffs[2,])
dev.off()
diffs<-do.call(rbind,lapply(fit$virusId,function(vir){
  do.call(rbind,lapply(fit$alleleId,function(xx){
    out<-data.frame(do.call(rbind,lapply(fit$alleleId,function(yy){
      if(xx==yy)return(c(0,0,.5))
      if(xx==1)return(crIMean(-mat[,sprintf('alleleVirus[%d,%d]',yy,vir)]))
      if(yy==1)return(crIMean(mat[,sprintf('alleleVirus[%d,%d]',xx,vir)]))
      return(crIMean(mat[,sprintf('alleleVirus[%d,%d]',xx,vir)]-mat[,sprintf('alleleVirus[%d,%d]',yy,vir)]))
    })))
    out$allele1<-xx
    out$allele2<-fit$alleleId
    out$virus<-vir
    return(out)
  }))
}))
diffs<-diffs[diffs$allele1<diffs$allele2,]
diffs$vir<-names(fit$virusId)[diffs$vir]
diffs$all1<-names(fit$alleleId)[diffs$allele1]
diffs$all2<-names(fit$alleleId)[diffs$allele2]


pdf('monkDiff.pdf',height=10,width=8)
sig01<-diffs[diffs$V3<.005|diffs$V3>.995,]
par(mar=c(5,7,.1,.6))
virOrder<-rev(unlist(lapply(virLookup,names)))
tmp<-means[!grepl('Mock',rownames(means)),]
tmp2<-tmp[virOrder,alleleLookup]
rownames(tmp2)<-do.call(c,unname(virLookup))[rownames(tmp2)]
sig01$row<-structure(1:nrow(tmp2),.Names=names(rownames(tmp2)))[sig01$vir]
sig01$col1<-structure(1:ncol(tmp2),.Names=colnames(tmp2))[sig01$all1]
sig01$rowOffset<-ave(sig01$row,sig01$row,FUN=function(xx)xx+seq(-length(xx)/50,length(xx)/65,length.out=length(xx)))+ave(sig01$col1,sig01$row,FUN=function(xx)cumsum(c(FALSE,xx[-1]!=xx[-length(xx)])*.02))
sig01$col2<-structure(1:ncol(tmp2),.Names=colnames(tmp2))[sig01$all2]
plotFunc(tmp2)
abline(h=cumsum(rev(unlist(sapply(virLookup,length))))[-length(virLookup)]+.5,lwd=2)
segments(sig01$col1,sig01$rowOffset,sig01$col2,sig01$rowOffset,col='#00000099')
points(c(sig01$col1,sig01$col2),c(sig01$rowOffset,sig01$rowOffset),col='#00000099',pch=21,bg='#00000099',cex=.3)
dev.off()
out<-diffs[,c('vir','all1','all2','V3','X2.5.','X97.5.')]
colnames(out)<-c('Virus','Allele1','Allele2','p(fold change<1)','Lower95CrI','Upper95CrI')
out$Lower95CrI<-exp(out$Lower95CrI)
out$Upper95CrI<-exp(out$Upper95CrI)
write.csv(out,'monkeyDiffs.csv')

stanCode2<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nAllele;
    int<lower=0> nRep;
    int<lower=0> nGroup;
    int<lower=0> nExperiment;
    vector[nRep] rep;
    int<lower=0> allele[nRep];
    int<lower=0> virus[nRep];
    int<lower=0> virusGroup[nVirus];
    int<lower=0> experiment[nRep];
    real backMean;
    real<lower=0> backSd;
    real<lower=0> tNu;
  }
  parameters {
    matrix[nVirus,nAllele] alleleVirus;
    matrix[nVirus,nExperiment-1] virusExperiment;
    matrix[nAllele,nGroup] alleleGroup;
    real<lower=0> sdGroup[nGroup];
    real<lower=0> repSd;
    vector<lower=0>[nRep] background;
    real<lower=0> alleleSd;
    vector[nAllele] alleleMeans;
  }
  transformed parameters{
    vector[nRep] beta;
    for(ii in 1:nRep){
      beta[ii]=alleleVirus[virus[ii],allele[ii]];
      if(experiment[ii]!=1)beta[ii]=beta[ii]+virusExperiment[virus[ii],experiment[ii]-1];
    }
  }
  model {
    rep~normal(log(exp(beta)+background),repSd);
    for(ii in 1:nVirus)alleleVirus[ii,]~student_t(tNu,alleleGroup[,virusGroup[ii]],sdGroup[virusGroup[ii]]);
    for(ii in 1:nGroup)alleleGroup[,ii]~student_t(tNu,alleleMeans,alleleSd);
    for(ii in 2:nExperiment)virusExperiment[,ii-1]~normal(0,5);
    repSd~gamma(1,.1);
    sdGroup~gamma(1,.1);
    background~normal(backMean,backSd);
    alleleMeans~normal(-2,3);
    //virusMeans~normal(0,10);
    alleleSd~gamma(1,.1);
  }
'
modAllele2 <- stan_model(model_code = stanCode2)

fit<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$siv,xx$rep,xx$date,mocks,modAllele2,chains=50,nIter=10000))
mat<-as.matrix(fit$fit)
pdf('monkDiff3.pdf',height=10,width=8)
plotFit(fit)
dev.off()

out<-diffs[,c('vir','all1','all2','V3','X2.5.','X97.5.')]
colnames(out)<-c('Virus','Allele1','Allele2','p(fold change<1)','Lower95CrI','Upper95CrI')
out$Lower95CrI<-exp(out$Lower95CrI)
out$Upper95CrI<-exp(out$Upper95CrI)
write.csv(out,'monkeyDiffs2.csv')

fit2<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$siv,xx$rep,xx$date,mocks,modAllele2,chains=50,nIter=10000))
save(fit2,file='tmp2.Rdat')
fit3<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$species,xx$rep,xx$date,mocks,modAllele2,chains=50,nIter=10000,tNu=10))
fitGor<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&!grepl('Bonobo|N15T',gor$allele),],fitModel(xx$virus,xx$allele,xx$siv,xx$rep,xx$date,mocks,modAllele2,chains=20,nIter=10000))
fitCpz<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&grepl('Bonobo|Human',gor$allele),],fitModel(xx$virus,xx$allele,xx$siv,xx$rep,xx$date,mocks,modAllele2,chains=20,nIter=10000))

print(fitGor$fit,c('beta','background'),include=FALSE)
print(fitCpz$fit,c('beta','background'),include=FALSE)

pdf('alleleGroup.pdf',width=14,height=5)
  par(mar=c(0,0,0,0))
  nAllele<-max(fit$alleleId)
  nGroup<-max(fit$groupId)
  layout(cbind(0,rbind(matrix(1:(nAllele*nGroup),nrow=nGroup),0)),height=c(rep(1,nGroup),.5),width=c(.5,rep(1,nAllele)))
  allGroup<-mat[,grepl('alleleGroup',colnames(mat))]
  allSd<-mat[,grepl('sdGroup',colnames(mat))]
  plusSd<-lapply(structure(colnames(allGroup),.Names=colnames(allGroup)),function(xx)density(LaplacesDemon::rst(nrow(allGroup),allGroup[,xx],allSd[,sub('alleleGroup\\[[0-9]+,','sdGroup[',xx)],3)))
  #xRange<-quantile(allGroup,c(.005,.995))
  xRange<-log(c(1e-6,.5))
  denses<-apply(allGroup,2,density)
  yLim<-c(0,max(sapply(denses,function(xx)max(xx$y))))
  for(jj in fit$alleleId){
    for(ii in fit$groupId){
      thisCol<-sprintf('alleleGroup[%d,%d]',jj,ii)
      thisDat<-denses[[thisCol]]
      thisVirusIds<-fit$virusId[names(fit$virusGroup)[fit$virusGroup==ii]]
      thisVirus<-mat[,sprintf('alleleVirus[%d,%d]',thisVirusIds,jj)]
      virDense<-apply(thisVirus,2,density)
      plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=exp(xRange),ylim=yLim,log='x',yaxs='i')
      #plot(1,1,type='n',xaxt='n',yaxt='n',xlab='',ylab='',xlim=exp(xRange),ylim=c(-yLim[2],yLim[2]),log='x')
      #polygon(exp(thisDat$x),thisDat$y,col='#0000FF33')
      #lapply(virDense,function(xx)polygon(exp(xx$x),-xx$y,col='#00000022',border='#00000011'))
      lapply(virDense,function(xx)polygon(exp(xx$x),xx$y,col='#00000011',border='#00000022'))
      polygon(exp(plusSd[[thisCol]]$x),plusSd[[thisCol]]$y,col='#FF000099')
      if(jj==1)mtext(names(fit$groupId)[ii],2,xpd=NA,line=.5)
      if(jj==1&ii==2)mtext('Probability density',2,xpd=NA,cex=1.2,line=3)
      if(ii==1)mtext(names(fit$alleleId)[jj],3,xpd=NA,line=-2)
      if(ii==nGroup)dnar::logAxis(1,axisVals=c(-6,-4,-2))
      if(ii==nGroup&jj==4)mtext('Proportion positive cells',1,xpd=NA,line=3)
    }
  }
dev.off()

save(fit,file='monkeyDiff.Rdat')


stanCode3<-'
  data {
    int<lower=0> nVirus;
    int<lower=0> nAllele;
    int<lower=0> nRep;
    int<lower=0> nGroup;
    int<lower=0> nExperiment;
    vector[nRep] rep;
    int<lower=0> allele[nRep];
    int<lower=0> virus[nRep];
    int<lower=0> virusGroup[nVirus];
    int<lower=0> experiment[nRep];
    real backMean;
    real<lower=0> backSd;
    real<lower=0> tNu;
  }
  parameters {
    matrix[nVirus,nAllele] alleleVirus;
    matrix[nVirus,nExperiment-1] virusExperiment;
    matrix[nAllele,nGroup] alleleGroup;
    real<lower=0> sdWithinSpecies;
    real<lower=0> repSd;
    vector<lower=0>[nRep] background;
    real<lower=0> alleleSd;
    vector[nAllele] sdMult;
    vector[nAllele] alleleMeans;
    vector[nVirus] virusMeans;
    //real<lower=0> alleleSdScale;
    //real<lower=0> alleleMeansSd;
    //real alleleMeansMean;
  }
  transformed parameters{
    vector[nRep] beta;
    for(ii in 1:nRep){
      beta[ii]=alleleVirus[virus[ii],allele[ii]]+virusMeans[virus[ii]];
      if(experiment[ii]!=1)beta[ii]=beta[ii]+virusExperiment[virus[ii],experiment[ii]-1];
      beta[ii]=logit(fmin(inv_logit(beta[ii])+background[ii],.9999));
    }
  }
  model {
    rep~normal(beta,repSd);
    //for(ii in 1:nVirus)alleleVirus[ii,]~student_t(tNu,alleleGroup[,virusGroup[ii]],sdWithinSpecies);
    for(ii in 1:nVirus)alleleVirus[ii,]~normal(alleleGroup[,virusGroup[ii]],sdWithinSpecies);
    for(ii in 1:nAllele)alleleGroup[ii,]~normal(alleleMeans[ii],alleleSd*exp(sdMult[ii]));
    for(ii in 2:nExperiment)virusExperiment[,ii-1]~normal(0,5);
    repSd~gamma(1,.1);
    alleleSd~gamma(1,.1);
    sdMult~normal(0,.693); //log(2)
    background~normal(backMean,backSd);
    alleleMeans~normal(-2,2.5);
    //alleleMeansMean~normal(-2,5);
    virusMeans~normal(0,2);
    sdWithinSpecies~gamma(1,.1);
    //alleleMeansSd~gamma(1,.1);
  }
'
modAllele3 <- stan_model(model_code = stanCode3)

fitMonk<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$species,xx$repNo0,xx$date,mocks,modAllele3,chains=50,nIter=10000,logFunc=logit))
print(fit6$fit,'alleleMeansSd')
pdf('monkDiff3.pdf',height=10,width=8)
plotFit(fit5)
dev.off()

fitGor<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&!grepl('Bonobo|N15T',gor$allele),],fitModel(xx$virus,xx$allele,xx$species,xx$rep,xx$date,mocks,modAllele3,chains=50,nIter=30000,logFunc=logit))
fitCpz<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&grepl('Bonobo|Human',gor$allele),],fitModel(xx$virus,xx$allele,xx$species,xx$rep,xx$date,mocks,modAllele3,chains=50,nIter=30000,logFunc=logit))

pdf('alleleGroup3.pdf',width=10,height=10)
plotAlleles(fitMonk,xRange=log(c(1e-9,1)),isT=FALSE,expFunc=boot::inv.logit,sdMult=TRUE)
dev.off()
mat<-as.matrix(fit5$fit)
alMeans<-mat[,grepl('alleleMeans\\[',colnames(mat))]
apply(alMeans,2,mean)
meanDiffs<-do.call(rbind,lapply(structure(colnames(alMeans),.Names=colnames(alMeans)),function(xx){do.call(rbind,lapply(structure(colnames(alMeans),.Names=colnames(alMeans)),function(yy){
  tmp<-crIMean(mat[,xx]-mat[,yy])
  pr<-mean((mat[,xx]-mat[,yy])<0)
  out<-data.frame('allele1'=xx,'allele2'=yy,'foldDiff'=exp(tmp[3]),'lowerCrI'=exp(tmp[1]),'upperCrI'=exp(tmp[2]),'probLess'=pr)
}))}))
meanDiffs$all1<-names(fit$alleleId)[as.numeric(sub('alleleMeans\\[([0-9]+)\\]','\\1',meanDiffs$allele1))]
meanDiffs$all2<-names(fit$alleleId)[as.numeric(sub('alleleMeans\\[([0-9]+)\\]','\\1',meanDiffs$allele2))]
rownames(meanDiffs)<-NULL
meanDiffs$diff<-apply(meanDiffs[,c('lowerCrI','upperCrI')],1,dnar::conservativeBoundary,base=1)

pdf('monkDiff4.pdf',height=10,width=8)
  plotFit(fitMonk,means)
dev.off()


zz<-rnorm(100000,mean(mocks/100),sd(mocks/100))
pdf('test.pdf');hist(mat[,'background[1]']);hist(mat[,'background[123]']);hist(mocks/100);hist(zz[zz>0]);dev.off()
