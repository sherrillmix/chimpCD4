source('functions.R')

mb<-read.csv('MB897 glycan mutants for Scott.csv',stringsAsFactors=FALSE)
colnames(mb)[1]<-'env'
colnames(mb)[2:ncol(mb)]<-sprintf('%s %d',rep(colnames(mb)[seq(2,ncol(mb),3)],each=3),1:3)
rownames(mb)<-mb$env

pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'))

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair)compareAlleles(mod,pair[1],pair[2],mb))
  names(selectFits)<-sapply(pairs,paste,collapse='-')
}

xlims<-exp(range(unlist(lapply(selectFits,findLims))))
dummy<-lapply(names(selectFits),function(xx){message(xx);print(selectFits[[xx]],pars=c('metaFoldChange','foldChange','repEffect'))})
pdf('mbFits.pdf',width=4,height=4);par(mar=c(3.5,11,2,.4));lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(mb),main=xx,cols=envCols,xlims=xlims));dev.off()
pdf('mbRaw.pdf',width=8,height=6);par(mar=c(3,4,1,.2));lapply(pairs,function(xx)plotRaw(mb[,grep(xx[1],colnames(mb))],mb[,grep(xx[2],colnames(mb))],cols=envCols));dev.off()


