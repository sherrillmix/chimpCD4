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


meansGori<-tapply(gor$norm,gor[colnames(gor)!='Human',c('virus','allele')],mean,na.rm=TRUE)

virLookup<-list(
  'Other'=c('SIVmusGAB11'='SIVmusGAB11','SIVmus1-54'='SIVmus1-54','SIVgsn CN71'='SIVgsnCN71','SIVmon'='SIVmonCM99','SIVdeb'='SIVdeb04CMP','SIVlst7'='SIVlstUS7','SIVagmTAN1'='SIVagmTAN1','SIVwrcGM05'='SIVwrcGM05','SIVasc'='SIVasc','SIVsmmFTq'='SIVsmmFTq','SIVsmmSL92b'='SIVsmmSL92b'),
  'SIVcpz'=c('MT145'='SIVcpzMT145','MB897'='SIVcpzMB897','BF1167'='SIVcpzBF1167','TAN2'='SIVcpzTAN2'),
  'SIVgor'=c('SIVgorCP2135'='SIVgorCP2135','SIVgorBPID1'='SIVgorBPID1','SIVgorBQID2'='SIVgorBQID2')
)
alleleLookup<-unique(monk$allele)
means<-tapply(monk$rep,monk[,c('virus','allele')],mean,na.rm=TRUE)
mocks<-c(gor[grepl('Mock',gor$virus)&!is.na(gor$rep),'rep'],monk[monk$virus=='Mock'&!is.na(monk$rep),'rep'])


fitMonk<-dnar::withAs(xx=monk[monk$virus!='Mock'&!is.na(monk$rep),],fitModel(xx$virus,xx$allele,xx$species,xx$repNo0,xx$date,mocks,modAllele,chains=50,nIter=50000,thin=4,logFunc=logit))
print(fitMonk$fit,'repSd')
#fitGor<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&!grepl('Bonobo|N15T',gor$allele),],fitModel(xx$virus,xx$allele,xx$species,xx$rep,xx$date,mocks,modAllele,chains=50,nIter=50000,thin=4,logFunc=logit))
#print(fitGor$fit,'repSd')
fitGor<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&!grepl('Bonobo',gor$allele),],fitModel(xx$virus,ifelse(xx$allele=='AHSM N15T','AHSM',xx$allele),xx$species,xx$rep,xx$date,mocks,modAlleleWithMods,ifelse(xx$allele=='AHSM N15T','N15T','Unmodified'),chains=50,nIter=50000,thin=4,logFunc=logit))
print(fitGor$fit,'repSd')
fitCpz<-dnar::withAs(xx=gor[!grepl('Mock',gor$virus)&!is.na(gor$rep)&grepl('Bonobo|Human',gor$allele),],fitModel(xx$virus,xx$allele,xx$species,xx$rep,xx$date,mocks,modAllele,chains=50,nIter=50000,thin=4,logFunc=logit))
print(fitCpz$fit,'repSd')
if(!file.exists('work/fitMonk.Rdat'){
  message('Saving data')
  save(fitMonk,file='work/fitMonk.Rdat')
  save(fitGor,file='work/fitGor.Rdat')
  save(fitCpz,file='work/fitCpz.Rdat')
}else{
  message('Not overwriting save')
}

pdf('monkDiff3.pdf',height=10,width=8)
plotFit(fitMonk)
dev.off()

pdf('alleleGroup3.pdf',width=10,height=10)
plotAlleles(fitMonk,xRange=log(c(1e-9,1)),isT=FALSE,expFunc=boot::inv.logit,sdMult=TRUE)
dev.off()


getDiffs<-function(fit){
  mat<-as.matrix(fit$fit)
  diffs<-do.call(rbind,parallel::mclapply(c('General'=-99,fit$virusId),function(vir){
    do.call(rbind,lapply(fit$alleleId,function(xx){
      out<-data.frame(do.call(rbind,lapply(fit$alleleId,function(yy){
        if(xx==yy)stats<-c(0,0,0,.5)
        else if(vir==-99)stats<-crIMean(mat[,sprintf('alleleMeans[%d]',xx)]-mat[,sprintf('alleleMeans[%d]',yy)])
        else stats<-crIMean(mat[,sprintf('alleleVirus[%d,%d]',vir,xx)]-mat[,sprintf('alleleVirus[%d,%d]',vir,yy)])
        #if(out[3]<0){
          #tmp<-xx;xx<-yy;yy<-tmp
          #stats<-crIMean(mat[,sprintf('alleleVirus[%d,%d]',vir,xx)]-mat[,sprintf('alleleVirus[%d,%d]',vir,yy)])
        #}
        out<-data.frame('allele1'=xx,'allele2'=yy,'mean'=stats[3],'lower'=stats[1],'upper'=stats[2],'p'=stats[4],stringsAsFactors=FALSE)
        return(out)
      })))
      out$vir<-vir
      return(out)
    }))
  },mc.cores=20))
  diffs<-diffs[diffs$mean>0&diffs$allele1!=diffs$allele2,]
  if(nrow(diffs)!=choose(max(fit$alleleId),2)*(1+max(fit$virusId)))stop('Missing rows')
  diffs$virus[diffs$vir!=-99]<-names(fit$virusId)[diffs$vir[diffs$vir!=-99]]
  diffs$virus[diffs$vir==-99]<-'General'
  diffs$all1<-names(fit$alleleId)[diffs$allele1]
  diffs$all2<-names(fit$alleleId)[diffs$allele2]
  return(diffs)
}
monkDiffs<-getDiffs(fitMonk)
monkDiffs$group<-'Other'
gorDiffs<-getDiffs(fitGor)
gorDiffs$group<-'Gorilla'
bonoDiffs<-getDiffs(fitCpz)
bonoDiffs$group<-'Bonobo'
out<-outPrep<-rbind(monkDiffs,gorDiffs,bonoDiffs)[,c('group','virus','all1','all2','p','mean','lower','upper')]
rownames(out)<-NULL
#censor high values >1000 to lower CrI
outPrep$censor<-out$mean>log(1000)|((out$upper-out$lower)>5&out$lower>0)
out[outPrep$censor,'upper']<-out[outPrep$censor,'lower']
out[outPrep$censor,'mean']<-out[outPrep$censor,'lower']
sigFig<-function(xx,digits=2)sub('\\.$','',formatC(signif(xx,digits=digits), digits=digits,format='fg',flag="#"))
out$upper<-sigFig(exp(out$upper))
out$lower<-sigFig(exp(out$lower))
out$mean<-sigFig(exp(out$mean))
out$upper[outPrep$upper>log(1000)]<-'>1000*'
out$upper[out$mean==out$lower&outPrep$censor]<-sprintf('>%s*',(out$lower[out$mean==out$lower&outPrep$censor]))
out$mean[out$mean==out$lower&outPrep$censor]<-sprintf('>%s*',out$lower[out$mean==out$lower&outPrep$censor])
out$p<-sigFig(out$p)
out[outPrep$p<.00001,'p']<-'<0.00001'
colnames(out)<-c('Group','Virus','AlleleA','AlleleB','p(fold change < 1)','Estimated fold change','Lower 95% CrI','Upper 95% CrI')
write.csv(out,'out/allDiffs.csv',row.names=FALSE)
for(ii in unique(out$Group))write.csv(out[out$Group==ii,-1],sprintf('out/%sDiffs.csv',ii),row.names=FALSE)

mat<-as.matrix(fitGor$fit)
stats<-as.data.frame(t(apply(mat[,grep('[mM]od',colnames(mat))],2,crIMean)))
colnames(stats)<-c('lower','upper','estimate','p')
stats<-stats[rownames(stats)!='modSd[1]',]
stats$virus<-ifelse(rownames(stats)=='modMeans[1]','Overall',names(fitGor$virusId)[suppressWarnings(as.numeric(sub('virusMod\\[([0-9]+),.*','\\1',rownames(stats))))])
out<-stats[,c('virus','lower','upper','estimate','p')]
out[,c('lower','upper','estimate')]<-apply(exp(out[,c('lower','upper','estimate')]),2,sigFig)
out[,c('p')]<-apply(out[,c('p'),drop=FALSE],2,sigFig)
out[stats$p<.00001,'p']<-'<0.00001'
out$AlleleA<-'AHSM'
out$AlleleB<-'AHSM N15T'
out<-out[,c('virus','AlleleA','AlleleB','p','estimate','lower','upper')]
colnames(out)<-c('Virus','AlleleA','AlleleB','p(fold change < 1)','Estimated fold change','Lower 95% CrI','Upper 95% CrI')
write.csv(out,'out/N15T.csv',row.names=FALSE)

mat<-as.matrix(fitCpz$fit)
stats<-as.data.frame(rbind('Overall'=crIMean(mat[,'alleleMeans[3]']-mat[,'alleleMeans[2]']),do.call(rbind,lapply(fitCpz$virusId,function(xx)crIMean(mat[,sprintf('alleleVirus[%d,3]',xx)]-mat[,sprintf('alleleVirus[%d,2]',xx)])))))
colnames(stats)<-c('lower','upper','estimate','p')
stats$virus<-rownames(stats)
out<-stats[,c('virus','lower','upper','estimate','p')]
out[,c('lower','upper','estimate')]<-apply(exp(out[,c('lower','upper','estimate')]),2,sigFig)
out[,c('p')]<-apply(out[,c('p'),drop=FALSE],2,sigFig)
out[stats$p<.00001,'p']<-'<0.00001'
write.csv(out,'out/I83T83.csv',row.names=FALSE)
print(out,row.names=FALSE)

mat<-as.matrix(fitCpz$fit)
stats<-do.call(rbind,lapply(2:3,function(yy){
  out<-as.data.frame(rbind('Overall'=crIMean(mat[,'alleleMeans[1]']-mat[,sprintf('alleleMeans[%d]',yy)]),do.call(rbind,lapply(fitCpz$virusId,function(xx)crIMean(mat[,sprintf('alleleVirus[%d,1]',xx)]-mat[,sprintf('alleleVirus[%d,%d]',xx,yy)])))))
  out$allele<-names(fitCpz$alleleId)[yy]
  out$virus<-rownames(out)
  out
}))
colnames(stats)<-c('lower','upper','estimate','p','allele','virus')
out<-stats[,c('allele','virus','lower','upper','estimate','p')]
out[,c('lower','upper','estimate')]<-apply(exp(out[,c('lower','upper','estimate')]),2,sigFig)
out[,c('p')]<-apply(out[,c('p'),drop=FALSE],2,sigFig)
out[stats$p<.00001,'p']<-'<0.00001'
write.csv(out,'out/Human_vs_I83T83.csv',row.names=FALSE)
print(out,row.names=FALSE)
