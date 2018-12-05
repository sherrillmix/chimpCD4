source('functions.R')

back<-as.numeric(as.vector(read.csv('Background per allele for Scott.csv',stringsAsFactors=FALSE)[1,-1]))
#summary(glm(back~as.factor(rep(1:10,each=3))))
pdf('test.pdf');plot(1:length(back),back,xlab='Background sample',ylab='Background % infection');abline(v=seq(3,27,3)+.5);dev.off()

cpz<-read.csv('SIVcpz Envs for Scott.csv',stringsAsFactors=FALSE)
colnames(cpz)[1]<-'env'
colnames(cpz)[2:ncol(cpz)]<-sprintf('%s %d',rep(colnames(cpz)[seq(2,ncol(cpz),3)],each=3),1:3)
cpz[cpz$env=='Gab2','env']<-'GAB2'
rownames(cpz)<-cpz$env
cpz<-cpz[dnar::orderIn(rownames(cpz),names(envCols)),]


pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('Human','QQKVP'),c('QQKVP','QQKVT'),c('QQNVP','RQNVP'),c('QQNVP','QRNVP'),c('QQNIP','QRNIP'),c('QQKVT','QRKVT'),c('QQNVP','QQKVP'),c('QQNVT','QQKVT'),c('QQNVP','QQNIP'),c('QRNVP','QRNIP'))

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair)compareAlleles(mod,pair[1],pair[2],cpz))
  names(selectFits)<-sapply(pairs,paste,collapse='-')
}

xlims<-exp(range(unlist(lapply(selectFits,findLims))))
dummy<-lapply(names(selectFits),function(xx){message(xx);print(selectFits[[xx]],pars=c('metaFoldChange','foldChange','repEffect'))})
pdf('cpzFits.pdf',width=4,height=4);par(mar=c(3.5,6,2,.4));lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(cpz),main=xx,cols=envCols,xlims=xlims));dev.off()
pdf('cpzRaw.pdf',width=8,height=6);par(mar=c(3,4,1,.2));lapply(pairs,function(xx)plotRaw(cpz[,grep(xx[1],colnames(cpz))],cpz[,grep(xx[2],colnames(cpz))],cols=envCols));dev.off()


