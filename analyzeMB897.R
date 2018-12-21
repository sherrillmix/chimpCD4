set.seed(12347)
source('functions.R')

mb<-read.csv('MB897.csv',stringsAsFactors=FALSE)
colnames(mb)[1]<-'env'
colnames(mb)[2:ncol(mb)]<-sprintf('%s %d',rep(colnames(mb)[seq(2,ncol(mb),3)],each=3),1:3)
rownames(mb)<-mb$env

pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('Human','QQNVT'))
#from analyzeCpz.R
wtPar<-list('Human-QQNVP'=c(-.632,.365),'QQNVP-QQNVT'=c(-1.05,.259),'Human-QQNVT'=c(-1.68,.573))

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair){
    wtMuSd<-wtPar[[paste(pair,collapse='-')]]
    compareAlleles(modMB,pair[1],pair[2],mb,wtPar=wtMuSd)
  })
  names(selectFits)<-sapply(pairs,paste,collapse='-')
}

xlims<-exp(range(unlist(lapply(selectFits,findLims))))
pdf('out/mbFits.pdf',width=4,height=4)
  par(mar=c(3.5,11,2,.4))
  lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(mb),main=sub('-','/',xx),xlims=xlims,special='WT',xlab='Fold change from WT',cols=c("WT"='#BBBBBB')))
dev.off()

pdf('out/mbRaw.pdf',width=8,height=6)
  par(mar=c(3,4,1,.2))
  lapply(pairs,function(xx)plotRaw(mb[,grep(xx[1],colnames(mb))],mb[,grep(xx[2],colnames(mb))],cols=envCols))
dev.off()

stats<-lapply(names(selectFits),function(xx){
  out<-assignNames(pullRanges(selectFits[[xx]],convertFunc=exp),rownames(mb))
  data.frame('comparison'=xx,'stat'=rownames(out),out)
})
names(stats)<-names(selectFits)
write.csv(do.call(rbind,stats),'out/mbStats.csv')


