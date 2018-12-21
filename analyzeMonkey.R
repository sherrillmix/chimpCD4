source('functions.R')

monkey1<-read.csv('monkeyGlycan.csv',stringsAsFactors=FALSE)
monkey2<-read.csv('monkeyOther.csv',stringsAsFactors=FALSE)

colnames(monkey2)[1]<-colnames(monkey1)[1]<-'env'
colnames(monkey1)[2:ncol(monkey1)]<-sprintf('%s %d',rep(colnames(monkey1)[seq(2,ncol(monkey1),3)],each=3),1:3)
colnames(monkey2)[2:ncol(monkey2)]<-sprintf('%s %d',rep(colnames(monkey2)[seq(2,ncol(monkey2),3)],each=3),1:3)
rownames(monkey1)<-monkey1$env
rownames(monkey2)<-monkey2$env

monkey1<-monkey1[dnar::orderIn(rownames(monkey1),names(envCols)),]
monkey2<-monkey2[dnar::orderIn(rownames(monkey2),names(envCols)),]
envIn2<-rownames(monkey1) %in% rownames(monkey2)
if(any(rownames(monkey1)[envIn2]!=rownames(monkey2)))stop('Env order mixup')

pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('QQNVP','QQNVP.N32Q'),c('QQNVT','QQNVT.N32Q'),c('QQNVP','RQNVP'),c('QQNVP','QRNVP'),c('QQNVP','QQKVP'),c('QQNVP','QQNIP'))

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair){
    in1<-c(any(grepl(sprintf('%s [0-9]+',pair[1]),colnames(monkey1))),any(grepl(sprintf('%s [0-9]+',pair[1]),colnames(monkey2))))
    in2<-c(any(grepl(sprintf('%s [0-9]+',pair[2]),colnames(monkey1))),any(grepl(sprintf('%s [0-9]+',pair[2]),colnames(monkey2))))
    mixDat<-in1[1]!=in2[1]
    if(!mixDat)compareAlleles(modCpz,pair[1],pair[2],if(in1[1])monkey1 else monkey2)
    else compareAlleles(modMonkey,pair[1],pair[2],if(in1[1])monkey1[envIn2,] else monkey2,if(in2[1])monkey1[envIn2,] else monkey2)
  })
  names(selectFits)<-sapply(pairs,paste,collapse='-')
}


xlims2<-exp(range(unlist(lapply(selectFits,findLims))))
pdf('out/monkeyFits2.pdf',width=4,height=4)
par(mar=c(3.5,8,2,.4))
lapply(names(selectFits),function(xx){
  if(sum(grepl('foldChange',colnames(as.matrix(selectFits[[xx]]))))==9)envNames<-rownames(monkey1)
  else envNames<-rownames(monkey2)
  plotFit(selectFits[[xx]],envNames,main=sub('-','/',xx),cols=NULL,xlims=xlims2)
})
dev.off()

pdf('out/monkeyRaw.pdf',width=10,height=6)
par(mar=c(3,4,1,.2))
lapply(pairs,function(xx){
  message(xx[1],'-',xx[2])
  regex1<-sprintf('%s [0-9]+',xx[1])
  regex2<-sprintf('%s [0-9]+',xx[2])
  in1<-c(any(grepl(regex1,colnames(monkey1))),any(grepl(regex1,colnames(monkey2))))
  in2<-c(any(grepl(regex2,colnames(monkey1))),any(grepl(regex2,colnames(monkey2))))
  mixDat<-in1[1]!=in2[1]
  dat1<-if(in1[1]){if(mixDat)monkey1[envIn2,] else monkey1}else monkey2
  dat2<-if(in2[1]){if(mixDat)monkey1[envIn2,] else monkey1}else monkey2
  plotRaw(dat1[,grepl(regex1,colnames(dat1))],dat2[,grepl(regex2,colnames(dat2))],cols=envCols)
})
dev.off()

stats<-mclapply(names(selectFits),function(xx){
  if(sum(grepl('foldChange',colnames(as.matrix(selectFits[[xx]]))))==9)envNames<-rownames(monkey1)
  else envNames<-rownames(monkey2)
  out<-assignNames(pullRanges(selectFits[[xx]],convertFunc=exp),envNames)
  data.frame('comparison'=xx,'stat'=rownames(out),out)
})
nCols<-sapply(stats,ncol)
for(ii in unique(nCols)){
  out<-do.call(rbind,stats[nCols==ii])
  write.csv(out,sprintf('out/monkeyStats_%d.csv',ii-2),row.names=FALSE)
}

