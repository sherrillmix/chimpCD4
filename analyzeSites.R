
sites<-read.csv('All sites genotypes.csv',stringsAsFactors=FALSE)
sites<-sites[sites$Sample.ID!='GM4801',]
sites$siv<-ifelse(grepl('\\+',sites$SIV.infection),TRUE,ifelse(grepl('-',sites$SIV.infection),FALSE,NA))
sites$isLBMB<-sites$Site.code %in% c('LB','MB')
cols<-list('nt51'=colnames(sites)[grep('nt51',colnames(sites))],'25_40'=colnames(sites)[grep('25.and.40',colnames(sites))],'52_55_68'=colnames(sites)[grep('52.55.68',colnames(sites))])
sites$combo2A<-sprintf('%s%s',tolower(sites$exon2.nt51.allele.A),sites$exon2.amino.acids.25.and.40.allele.A)
sites$combo2B<-sprintf('%s%s',tolower(sites$exon2.nt51.allele.B),sites$exon2.amino.acids.25.and.40.allele.B)
cols2<-list('25_40'=colnames(sites)[grep('combo2',colnames(sites))],'52_55_68'=colnames(sites)[grep('52.55.68',colnames(sites))])
dummies<-do.call(cbind,lapply(names(cols),function(xx){
  pos<-strsplit(xx,'_')[[1]]
  allAllele<-unique(unlist(sites[,cols[[xx]]]))
  nChar<-nchar(allAllele[1])
  if(any(nchar(allAllele)!=nChar))stop('Different length alleles')
  out<-do.call(cbind,lapply(1:nChar,function(ii){
    chars<-apply(sites[,cols[[xx]]],1,function(yy)paste(substring(yy,ii,ii),collapse=''))
    uniqChar<-unique(unlist(strsplit(chars,'')))
    out<-do.call(cbind,lapply(uniqChar,function(yy)grepl(yy,chars)))
    colnames(out)<-sprintf('%d_%s',ii,uniqChar)
    return(out)
  }))
  for(ii in 1:length(pos))colnames(out)<-sub(sprintf('^%d_',ii),pos[ii],colnames(out))
  return(out)
}))
dummies2<-do.call(cbind,lapply(names(cols2),function(xx){
  allAllele<-unique(unlist(sites[,cols2[[xx]]]))
  out<-do.call(cbind,lapply(allAllele,function(ii){
      apply(sites[,cols2[[xx]]]==ii,1,any)
  }))
  colnames(out)<-sprintf('%s_%s',xx,allAllele)
  return(out)
}))

selector<-sites$isLBMB
summary(glm(sites$siv[selector]~.,data=as.data.frame(dummies[selector,]),family='binomial'))
selector<-sites$Site.code %in% c('LB','MB','GM') 
summary(glm(sites$siv[selector]~.+sites$isLBMB[selector],data=as.data.frame(dummies[selector,!grepl('nt',colnames(dummies))]),family='binomial'))
selector<-sites$Site.code %in% c('LB','MB','GM')
colSelect<-abs(apply(dummies2[selector,],2,mean)-.5)<.45#&!grepl('nt',colnames(dummies2))
mod<-as.data.frame(coef(summary(glm(sites$siv[selector]~.+sites$isLBMB[selector],data=as.data.frame(dummies2[selector,colSelect]),family='binomial'))))
colSelect<-abs(apply(dummies[selector,],2,mean)-.5)<.45#&!grepl('nt',colnames(dummies))
summary(glm(sites$siv[selector]~.+sites$isLBMB[selector],data=as.data.frame(dummies[selector,colSelect]),family='binomial'))


#Combined
selector<-sites$Site.code %in% c('LB','MB','GM')
colSelect<-apply(dummies2[selector,],2,sum)>0&apply(dummies2[selector,],2,sum)<sum(selector)
#colSelect<-colSelect & !colnames(dummies2) %in% c('25_40_QQ','52_55_68_KVP')
mod<-as.data.frame(coef(print(summary(glm(sites$siv[selector]~.+sites$isLBMB[selector],data=as.data.frame(dummies2[selector,colSelect]),family='binomial')))))
# LB/MB
selector<-sites$Site.code %in% c('LB','MB') 
colSelect<-apply(dummies[selector,],2,sum)>0&apply(dummies[selector,],2,sum)<sum(selector)
#avoid convergence problem with small baseline
#colSelect<-colSelect&colnames(dummies)!='nt51A'
mod_LM<-as.data.frame(coef(print(summary(glm(sites$siv[selector]~.,data=as.data.frame(dummies[selector,colSelect]),family='binomial')))))
# Gombe
selector<-sites$Site.code %in% c('GM') 
colSelect<-apply(dummies[selector,],2,sum)>0&apply(dummies[selector,],2,sum)<sum(selector)
#colSelect<-colSelect&colnames(dummies)!='25Q'
mod_G<-as.data.frame(coef(print(summary(glm(sites$siv[selector]~.,data=as.data.frame(dummies[selector,colSelect]),family='binomial')))))

#LB/MB allele levels
selector<-sites$Site.code %in% c('LB','MB')
colSelect<-apply(dummies2[selector,],2,sum)>0&apply(dummies2[selector,],2,sum)<sum(selector)
#colSelect<-colSelect & !colnames(dummies2) %in% c('nt51_A','52_55_68_KVP')
mod_LM_a<-as.data.frame(coef(print(summary(glm(sites$siv[selector]~.,data=as.data.frame(dummies2[selector,colSelect]),family='binomial')))))
#Gombe allele level
selector<-sites$Site.code %in% c('GM')
colSelect<-apply(dummies2[selector,],2,sum)>0&apply(dummies2[selector,],2,sum)<sum(selector)
#colSelect<-colSelect & !colnames(dummies2) %in% c('25_40_QQ')
mod_G_a<-as.data.frame(coef(print(summary(glm(sites$siv[selector]~.,data=as.data.frame(dummies2[selector,colSelect]),family='binomial')))))



cleanNames<-function(xx)sub('nt51-?G','nt51g',sub('nt51-?A','nt51a',gsub('_','-',sub('sites\\$isLBMB\\[selector\\]','inLBMB',gsub('`','',sub("TRUE",'',rownames(xx)))))))

rownames(mod)<-cleanNames(mod)
rownames(mod_LM)<-cleanNames(mod_LM)
rownames(mod_G)<-cleanNames(mod_G)
rownames(mod)<-sub('25-40-|52-55-68-','',rownames(mod))
rownames(mod_LM_a)<-cleanNames(mod_LM_a)
rownames(mod_G_a)<-cleanNames(mod_G_a)
rownames(mod_LM_a)<-sub('25-40-|52-55-68-','',rownames(mod_LM_a))
rownames(mod_G_a)<-sub('25-40-|52-55-68-','',rownames(mod_G_a))

#mod<-mod[rownames(mod)!='inLBMB',]
#mod_LM<-mod_LM[rownames(mod_LM)!='inLBMB',]
#mod_G<-mod_G[rownames(mod_G)!='inLBMB',]
mod<-mod[!rownames(mod) %in% c('(Intercept)','gRQ','KVP','inLBMB'),]
mod_LM<-mod_LM[!rownames(mod_LM) %in% c('(Intercept)','nt51a','inLBMB'),]
mod_G<-mod_G[!rownames(mod_G) %in% c('(Intercept)','25Q','inLBMB'),]
mod_LM_a<-mod_LM_a[!rownames(mod_LM_a) %in% c('(Intercept)','aQQ','KVP','inLBMB'),]
mod_G_a<-mod_G_a[!rownames(mod_G_a) %in% c('(Intercept)','gRQ','inLBMB'),]


plotGlm<-function(mod,xlim=range(c(mod$upper,mod$lower)),main='',xlab='Fold change in odds of SIVcpz infection',cex.axis=.9,...){
  mod$upper<-mod$Estimate+mod[,'Std. Error']*2
  mod$lower<-mod$Estimate-mod[,'Std. Error']*2
  ylims<-c(.5,nrow(mod))
  plot(1,1,type='n',xlim=exp(xlim),ylim=ylims,ylab='',xlab=xlab,yaxt='n',log='x',xaxt='n',mgp=c(2.25,1,0),bty='l',...)
  par(lheight=.7)
  mtext(sub(' ','\n',main),3,xpd=NA,at=ifelse(nchar(main)>20,1,1),line=ifelse(grepl(' ',main),0,.7),cex=cex.axis/1.2)
  points(exp(mod$Estimate),nrow(mod):1,pch=21,bg='black',cex=.8)
  segments(exp(mod$lower),nrow(mod):1,exp(mod$upper),nrow(mod):1)
  if(max(nchar(rownames(mod)))>100)dnar::slantAxis(2,nrow(mod):1,rownames(mod),location=.5,axisArgs=list(tcl=-.2),cex=cex.axis)
  else axis(2,nrow(mod):1,rownames(mod),las=1,xpd=NA,tcl=-.2,mgp=c(3,.5,0),cex.axis=cex.axis)
  dnar::logAxis(1,mgp=c(3,.4,0),cex.axis=cex.axis)
  abline(v=1,lty=2)
}
xlim<-dnar::withAs(tmp=rbind(mod,mod_LM,mod_G),range(c(tmp$Estimate+2*tmp[,'Std. Error'],tmp$Estimate-2*tmp[,'Std. Error'])))
pdf('glm.pdf',height=2.2,width=3.42)
par(mar=c(3.5,0.01,2.1,0.01),cex=.3)
layout(t(c(0,1,0,2,0,3,0)),width=c(2.4,5,2,5,3.8,5,.1))
plotGlm(mod_G[!grepl('Intercept',rownames(mod_G)),],main='Gombe (GM)',xlim=xlim,xlab='')
text(grconvertX(-.44, "nfc", "user"),grconvertY(.99, "nfc", "user"),'A',xpd=NA,cex=2,adj=0:1)
plotGlm(mod_LM[!grepl('Intercept',rownames(mod_LM)),],main='Lobéké/Manbélé (LB/MB)',xlim=xlim,xpd=NA)
text(grconvertX(-.48, "nfc", "user"),grconvertY(.99, "nfc", "user"),'B',xpd=NA,cex=2,adj=0:1)
plotGlm(mod[!grepl('Intercept',rownames(mod)),],main='Combined',xlim=xlim,xlab='')
text(grconvertX(-.48, "nfc", "user"),grconvertY(.99, "nfc", "user"),'C',xpd=NA,cex=2,adj=0:1)
exon2<-range(grep('[RQ][RQ]',rev(rownames(mod))))
exon3<-range(grep('[NK][VI][TP]',rev(rownames(mod))))
segments(grconvertX(-.5, "nfc", "user"),exon2[1]-.4,grconvertX(-.5, "nfc", "user"),exon2[2]+.4,xpd=NA)
segments(grconvertX(-.5, "nfc", "user"),exon2+c(-.4,.4),grconvertX(-.4, "nfc", "user"),exon2+c(-.4,.4),xpd=NA)
mtext('Exon 2',2,cex=.7,line=2.9,at=mean(exon2))
segments(grconvertX(-.5, "nfc", "user"),exon3[1]-.4,grconvertX(-.5, "nfc", "user"),exon3[2]+.4,xpd=NA)
segments(grconvertX(-.5, "nfc", "user"),exon3+c(-.4,.4),grconvertX(-.4, "nfc", "user"),exon3+c(-.4,.4),xpd=NA)
mtext('Exon 3',2,cex=.7,line=2.9,at=mean(exon3))
dev.off()

pdf('glm2.pdf',height=2,width=3.42)
par(mar=c(3.5,0.01,2.1,0.01))
layout(t(c(0,1,0,2,0,3,0,4,0)),width=c(2.7,5,3.5,5,2,5,3.8,5,1))
plotGlm(mod_G[!grepl('Intercept',rownames(mod_G)),],main='Gombe polymorphism',xlim=xlim,xlab='',cex.axis=.7)
#plotGlm(mod[!grepl('Intercept',rownames(mod)),],main='Combined',xlim=xlim,xlab='')
plotGlm(mod_G_a[!grepl('Intercept',rownames(mod_G_a)),],main='Gombe allele',xlim=xlim,xlab='',cex.axis=.7)
exon2<-range(grep('[RQ][RQ]',rev(rownames(mod_G_a))))
exon3<-range(grep('[NK][VI][TP]',rev(rownames(mod_G_a))))
segments(grconvertX(-.5, "nfc", "user"),exon2[1]-.4,grconvertX(-.5, "nfc", "user"),exon2[2]+.4,xpd=NA)
segments(grconvertX(-.5, "nfc", "user"),exon2+c(-.4,.4),grconvertX(-.4, "nfc", "user"),exon2+c(-.4,.4),xpd=NA)
mtext('Exon 2',2,cex=.7,line=2,at=mean(exon2))
segments(grconvertX(-.5, "nfc", "user"),exon3[1]-.4,grconvertX(-.5, "nfc", "user"),exon3[2]+.4,xpd=NA)
segments(grconvertX(-.5, "nfc", "user"),exon3+c(-.4,.4),grconvertX(-.4, "nfc", "user"),exon3+c(-.4,.4),xpd=NA)
mtext('Exon 3',2,cex=.7,line=2,at=mean(exon3))
plotGlm(mod_LM[!grepl('Intercept',rownames(mod_LM)),],main='Lobéké/Manbélé polymorphism',xlim=xlim,xpd=NA,cex.axis=.7)
plotGlm(mod_LM_a[!grepl('Intercept',rownames(mod_LM_a)),],main='Lobéké/Manbélé allele',xlim=xlim,xlab='',cex.axis=.7)
exon2<-range(grep('[RQ][RQ]',rev(rownames(mod_LM_a))))
exon3<-range(grep('[NK][VI][TP]',rev(rownames(mod_LM_a))))
segments(grconvertX(-.5, "nfc", "user"),exon2[1]-.4,grconvertX(-.5, "nfc", "user"),exon2[2]+.4,xpd=NA)
segments(grconvertX(-.5, "nfc", "user"),exon2+c(-.4,.4),grconvertX(-.4, "nfc", "user"),exon2+c(-.4,.4),xpd=NA)
mtext('Exon 2',2,cex=.7,line=2,at=mean(exon2))
segments(grconvertX(-.5, "nfc", "user"),exon3[1]-.4,grconvertX(-.5, "nfc", "user"),exon3[2]+.4,xpd=NA)
segments(grconvertX(-.5, "nfc", "user"),exon3+c(-.4,.4),grconvertX(-.4, "nfc", "user"),exon3+c(-.4,.4),xpd=NA)
mtext('Exon 3',2,cex=.7,line=2,at=mean(exon3))
dev.off()

pdf('glm3.pdf',height=2,width=3.42)
par(mar=c(3.5,2.1,2,0.2),mfrow=c(1,4))
plotGlm(mod_G[!grepl('Intercept',rownames(mod_G)),],main='Gombe polymorphism',xlim=xlim,xlab='',cex.axis=.7)
plotGlm(mod_G_a[!grepl('Intercept',rownames(mod_G_a)),],main='Gombe allele',xlim=xlim,xlab='',cex.axis=.7)
plotGlm(mod_LM[!grepl('Intercept',rownames(mod_LM)),],main='Lobéké/Manbélé polymorphism',xlim=xlim,xpd=NA,cex.axis=.7)
plotGlm(mod_LM_a[!grepl('Intercept',rownames(mod_LM_a)),],main='Lobéké/Manbélé allele',xlim=xlim,xlab='',cex.axis=.7)
dev.off()


write.csv(cbind(sites,dummies,dummies2),'sites.csv',row.names=FALSE)

selector<-sites$Site.code %in% c('LB','MB','GM')
do.call(cbind,lapply(structure(c('LB','MB','GM'),.Names=c('LB','MB','GM')),function(xx)apply(cbind(dummies2[sites$Site.code==xx,],'Total'=TRUE),2,sum)))
