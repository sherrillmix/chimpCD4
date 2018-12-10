source('functions.R')

mb<-read.csv('MB897 glycan mutants for Scott.csv',stringsAsFactors=FALSE)
colnames(mb)[1]<-'env'
colnames(mb)[2:ncol(mb)]<-sprintf('%s %d',rep(colnames(mb)[seq(2,ncol(mb),3)],each=3),1:3)
rownames(mb)<-mb$env

pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('Human','QQNVT'))
#from analyzeCpz.R
wtPar<-list('Human-QQNVP'=c(-.624,.358),'QQNVP-QQNVT'=c(-1.02,.223),'Human-QQNVT'=c(-1.64,.543))

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair){
    wtMuSd<-wtPar[[paste(pair,collapse='-')]]
    compareAlleles(mod3,pair[1],pair[2],mb,wtPar=wtMuSd)
  })
  names(selectFits)<-sapply(pairs,paste,collapse='-')
  selectFits2<-lapply(pairs,function(pair){
    wtMuSd<-wtPar[[paste(pair,collapse='-')]]
    compareAlleles(mod6,pair[1],pair[2],mb,wtPar=wtMuSd)
  })
  names(selectFits2)<-sapply(pairs,paste,collapse='-')
}

xlims<-exp(range(unlist(lapply(selectFits,findLims))))
dummy<-lapply(names(selectFits),function(xx){message(xx);print(selectFits[[xx]],pars=c('metaFoldChange','foldChange'))})
pdf('mbFits.pdf',width=4,height=4);par(mar=c(3.5,11,2,.4));lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(mb),main=sub('-','/',xx),cols=envCols,xlims=xlims,special='WT',xlab='Fold change from WT'));dev.off()
pdf('mbRaw.pdf',width=8,height=6);par(mar=c(3,4,1,.2));lapply(pairs,function(xx)plotRaw(mb[,grep(xx[1],colnames(mb))],mb[,grep(xx[2],colnames(mb))],cols=envCols));dev.off()


