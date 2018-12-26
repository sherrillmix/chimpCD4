set.seed(12345)
source('functions.R')

back<-as.numeric(as.vector(read.csv('data/background.csv',stringsAsFactors=FALSE)[1,-1]))
print(summary(log(back)))

cpz<-read.csv('data/SIVcpz.csv',stringsAsFactors=FALSE)
colnames(cpz)[1]<-'env'
colnames(cpz)[2:ncol(cpz)]<-sprintf('%s %d',rep(colnames(cpz)[seq(2,ncol(cpz),3)],each=3),1:3)
rownames(cpz)<-cpz$env

cpz<-cpz[order(sapply(rownames(cpz),function(xx)which(names(envCols)==xx))),]


pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('Human','QQKVP'),c('QQKVP','QQKVT'),c('QQNVP','RQNVP'),c('QQNVP','QRNVP'),c('QQNIP','QRNIP'),c('QQKVT','QRKVT'),c('QQNVP','QQKVP'),c('QQNVT','QQKVT'),c('QQNVP','QQNIP'),c('QRNVP','QRNIP'),c('Human','QQNVT'))

sameScale<-c('QQNVP-RQNVP','QQNVP-QRNVP','QQNIP-QRNIP','QQKVT-QRKVT')

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair)compareAlleles(modCpz,pair[1],pair[2],cpz))
  names(selectFits)<-sapply(pairs,paste,collapse='-')
}

pdf('out/cpzRaw.pdf',width=8,height=6)
  par(mar=c(3,4,1,.2))
  lapply(pairs,function(xx)plotRaw(cpz[,grep(xx[1],colnames(cpz))],cpz[,grep(xx[2],colnames(cpz))],cols=envCols))
dev.off()

xlims<-exp(range(unlist(lapply(selectFits[names(selectFits) %in% sameScale],findLims))))
xlims[1]<-max(xlims[1],1e-5)
xlims2<-exp(range(unlist(lapply(selectFits[!names(selectFits) %in% sameScale],findLims))))
pdf('out/cpzFits.pdf',width=4,height=4)
  par(mar=c(3.5,6,2,.4))
  lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(cpz),main=sub('-','/',xx),cols=NULL,xlims=if(xx %in% sameScale)xlims else xlims2))
dev.off()


apply(as.matrix(selectFits[['Human-QQNVP']])[,c('metaFoldChange','metaSd')],2,mean)
apply(as.matrix(selectFits[['QQNVP-QQNVT']])[,c('metaFoldChange','metaSd')],2,mean)
apply(as.matrix(selectFits[['Human-QQNVT']])[,c('metaFoldChange','metaSd')],2,mean)

stats<-lapply(names(selectFits),function(xx){
  out<-assignNames(pullRanges(selectFits[[xx]],convertFunc=exp),rownames(cpz))
  data.frame('comparison'=xx,'stat'=rownames(out),out)
})
names(stats)<-names(selectFits)
write.csv(do.call(rbind,stats),'out/cpzStats.csv')
