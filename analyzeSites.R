source('functions.R')

sites<-read.csv('data/siteGenotypes.csv',stringsAsFactors=FALSE)
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


ancestral<-c('25_40_gQQ','52_55_68_NVP','52_55_68_NVP','nt51G','25Q','40Q','52N','55V','68P','gQQ-gQQ','NVP-NVP')
# LB/MB
selector<-sites$Site.code %in% c('LB','MB') 
colSelect<-apply(dummies[selector,],2,sum)>0&apply(dummies[selector,],2,sum)<sum(selector)&!colnames(dummies) %in% ancestral
glm_LM<-glm(sites$siv[selector]~.,data=as.data.frame(dummies[selector,colSelect]),family='binomial')
print(anova(glm_LM,test='Chisq'))
mod_LM<-as.data.frame(coef(print(summary(glm_LM))))
# Gombe
selector<-sites$Site.code %in% c('GM') 
colSelect<-apply(dummies[selector,],2,sum)>0&apply(dummies[selector,],2,sum)<sum(selector)&!colnames(dummies) %in% ancestral
glm_G<-glm(sites$siv[selector]~.,data=as.data.frame(dummies[selector,colSelect]),family='binomial')
print(anova(glm_G,test='Chisq'))
mod_G<-as.data.frame(coef(print(summary(glm_G))))

#LB/MB allele levels
selector<-sites$Site.code %in% c('LB','MB')
colSelect<-apply(dummies2[selector,],2,sum)>0&apply(dummies2[selector,],2,sum)<sum(selector)&!colnames(dummies2) %in% ancestral
glm_LM_a<-glm(sites$siv[selector]~.,data=as.data.frame(dummies2[selector,colSelect]),family='binomial')
print(anova(glm_LM_a,test='Chisq'))
mod_LM_a<-as.data.frame(coef(print(summary(glm_LM_a))))
#Gombe allele level
selector<-sites$Site.code %in% c('GM')
colSelect<-apply(dummies2[selector,],2,sum)>0&apply(dummies2[selector,],2,sum)<sum(selector)&!colnames(dummies2) %in% ancestral
glm_G_a<-glm(sites$siv[selector]~.,data=as.data.frame(dummies2[selector,colSelect]),family='binomial')
print(anova(glm_G_a,test='Chisq'))
mod_G_a<-as.data.frame(coef(print(summary(glm_G_a))))

cleanNames<-function(xx)sub('nt51-?G','nt51g',sub('nt51-?A','nt51a',gsub('_','-',sub('sites\\$isLBMB\\[selector\\]','inLBMB',gsub('`','',sub("TRUE",'',rownames(xx)))))))

rownames(mod_LM)<-cleanNames(mod_LM)
rownames(mod_G)<-cleanNames(mod_G)
rownames(mod_LM_a)<-cleanNames(mod_LM_a)
rownames(mod_G_a)<-cleanNames(mod_G_a)
rownames(mod_LM_a)<-sub('25-40-|52-55-68-','',rownames(mod_LM_a))
rownames(mod_G_a)<-sub('25-40-|52-55-68-','',rownames(mod_G_a))

mod_LM<-mod_LM[!rownames(mod_LM) %in% c('(Intercept)','nt51a','inLBMB'),]
mod_G<-mod_G[!rownames(mod_G) %in% c('(Intercept)','25Q','inLBMB'),]
mod_LM_a<-mod_LM_a[!rownames(mod_LM_a) %in% c('(Intercept)','aQQ','KVP','inLBMB'),]
mod_G_a<-mod_G_a[!rownames(mod_G_a) %in% c('(Intercept)','gRQ','inLBMB'),]


xlim<-with(rbind(mod_LM,mod_G,mod_LM_a,mod_G_a),range(c(Estimate+2*`Std. Error`,Estimate-2*`Std. Error`)))

pdf('out/glm.pdf',height=2,width=3.42)
par(mar=c(3.5,0.01,2.1,0.01))
layout(t(c(0,1,0,2,0,3,0,4,0)),width=c(2.7,5,3.5,5,2,5,3.8,5,1))
plotGlm(mod_G[!grepl('Intercept',rownames(mod_G)),],main='Gombe polymorphism',xlim=xlim,xlab='',cex.axis=.7)
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


