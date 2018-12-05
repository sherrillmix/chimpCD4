envCd4<-read.csv('Monkey SIV Envs vs CD4 for Scott.csv')
colnames(envCd4)[1]<-'env'
colnames(envCd4)[2:ncol(envCd4)]<-sprintf('%s %d',rep(colnames(envCd4)[seq(2,ncol(envCd4),3)],each=3),1:3)
rownames(envCd4)<-envCd4$env
envCd4<-envCd4[,-1:-4]

pdf('heat.pdf',width=8,height=8)
  cols<-rev(heat.colors(100))
  breaks<-seq(-.00001,max(envCd4)+.0001,length.out=101)
  par(mar=c(5.5,5,.1,1))
  image(1:nrow(envCd4),1:ncol(envCd4),as.matrix(envCd4),xaxt='n',yaxt='n',xlab='',ylab='',breaks=breaks,col=cols)
  dnar::slantAxis(1,1:nrow(envCd4),rownames(envCd4))
  axis(2,1:ncol(envCd4),colnames(envCd4),las=1)
  box()
  dnar::insetScale(breaks,cols,main='Proportion of human')
  par(family='mono')
  heatmap(as.matrix(envCd4),col=cols,breaks=breaks,scale='none',margins=c(8,10))
  par(family='sans')
  dnar::insetScale(breaks,cols,main='Proportion of human')
  #heatmap(log(as.matrix(envCd4)),col=cols,scale='none',margins=c(8,10))
dev.off()

pdf('pca.pdf')
biplot(prcomp(envCd4),col=c('black','#FF000055'))
tsne<-Rtsne::Rtsne(envCd4[,-1:-4],perplexity=1.7)
tsnel<-Rtsne::Rtsne(log(envCd4[,-1:-4]),perplexity=1.7)
plot(tsne$Y,type='n')
text(tsne$Y,rownames(envCd4),xlab='t-SNE 1',ylab='t-SNE 2')
#plot(tsnel$Y,type='n')
#text(tsnel$Y,rownames(envCd4),xlab='t-SNE 1',ylab='t-SNE 2')
dev.off()

alleles<-unique(sapply(strsplit(colnames(envCd4),' '),'[[',1))
ps<-do.call(rbind,lapply(alleles[alleles!='QQNVP'],function(xx){
  sapply(rownames(envCd4),function(yy){
    #t.test(log(envCd4[yy,grep(xx,colnames(envCd4))]),mu=log(100))$p.value
    t.test(log(envCd4[yy,grep(xx,colnames(envCd4))]),log(envCd4[yy,grep('QQNVP',colnames(envCd4))]))$p.value
  })
}))
ps[,]<-p.adjust(ps,method='fdr')
rownames(ps)<-alleles[alleles!='QQNVP']
#ci<-do.call(rbind,lapply(alleles,function(xx){
  #sapply(rownames(envCd4),function(yy){
    #exp(dnar::conservativeBoundary(t.test(log(envCd4[yy,grep(xx,colnames(envCd4))]),mu=log(100))$conf.int,log(100)))
  #})
#}))

stack<-data.frame('env'=rep(rownames(envCd4),ncol(envCd4)),'alleleRep'=rep(colnames(envCd4),each=nrow(envCd4)),'infect'=unlist(envCd4))
stack$allele<-sub(' .*','',stack$alleleRep)
stack[,sprintf('aa%d',1:5)]<-do.call(rbind,strsplit(stack$allele,''))
base<-'QQNVP'
for(ii in 1:5){
  thisCol<-sprintf('aa%d',ii)
  baseAA<-substring(base,ii,ii)
  thisLevels<-c(baseAA,unique(stack[,thisCol][stack[,thisCol]!=baseAA]))
  stack[,sprintf('aa%d',ii)]<-factor(stack[,sprintf('aa%d',ii)],levels=thisLevels)
}


mod<-lm(log10(infect)~env+env*(aa1+aa2+aa3+aa4+aa5)+0,stack)
#mod<-lm(log10(infect)~aa1+aa2+aa3+aa4+aa5+0,stack[stack$env=='SIVmac251',])
xx<-as.data.frame(coef(summary(mod)))
xx$aa<-sapply(strsplit(rownames(xx),':'),function(yy){if(length(yy)>1)yy[2] else NA})
baseAa<-xx[grep('^aa[0-9][A-Z]$',rownames(xx)),'Estimate']
names(baseAa)<-rownames(xx)[grep('^aa[0-9][A-Z]$',rownames(xx))]
xx$Estimate<-xx$Estimate+ifelse(!is.na(xx$aa),baseAa[xx$aa],0)
sigs<-xx[xx[,'Pr(>|t|)']<.01&grepl('aa',rownames(xx)),]
rownames(sigs)[grep('^aa[0-9][A-Z]$',rownames(sigs))]<-sprintf('SIVagm:%s',rownames(sigs)[grep('^aa[0-9][A-Z]$',rownames(sigs))])
rownames(sigs)<-sub('envS','S',rownames(sigs))
#sigs[,'Estimate']<-exp(sigs[,'Estimate'])




pdf('ps.pdf',width=8,height=8)
  cols<-(heat.colors(100))
  breaks<-seq(min(log10(ps))-.0001,0,length.out=101)
  par(mar=c(5.5,7,.1,1))
  image(1:nrow(ps),1:ncol(ps),log10(as.matrix(ps)),xaxt='n',yaxt='n',xlab='',ylab='',breaks=breaks,col=cols)
  dnar::slantAxis(1,1:nrow(ps),rownames(ps))
  axis(2,1:ncol(ps),colnames(ps),las=1)
  box()
  dnar::insetScale(breaks,cols,at=c(-2,-1,0),label=c(.01,.1,1),main='p-value')
  points(which(ps<.05,arr.ind=TRUE),pch='*',col='#00000044')
dev.off()


