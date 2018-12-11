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


pairs<-list(c('Human','QQNVP'),c('QQNVP','QQNVT'),c('Human','QQKVP'),c('QQKVP','QQKVT'),c('QQNVP','RQNVP'),c('QQNVP','QRNVP'),c('QQNIP','QRNIP'),c('QQKVT','QRKVT'),c('QQNVP','QQKVP'),c('QQNVT','QQKVT'),c('QQNVP','QQNIP'),c('QRNVP','QRNIP'),c('Human','QQNVT'))

sameScale<-c('QQNVP-RQNVP','QQNVP-QRNVP','QQNIP-QRNIP','QQKVT-QRKVT')

if(!exists('selectFits')){
  selectFits<-lapply(pairs,function(pair)compareAlleles(mod,pair[1],pair[2],cpz))
  names(selectFits)<-sapply(pairs,paste,collapse='-')
  initF<-function(xx){
    out<-list(  
      metaFoldChange=0,
      metaSd=.1,
      prop1BaseRaw=matrix(rep(-2,3*nrow(cpz)),nrow=nrow(cpz)),
      metaRaw=rep(-1,nrow(cpz)),
      sigma=.1,
      repEffect=c(0,0),
      baseSigma=.1,
      baseMu=rep(-1,nrow(cpz)),
      back1=matrix(rnorm(3*nrow(cpz),-6.5,1),nrow=nrow(cpz)),
      back2=matrix(rnorm(3*nrow(cpz),-6.5,1),nrow=nrow(cpz)),
      trueProp1Raw=matrix(rnorm(3*nrow(cpz),0,1),nrow=nrow(cpz)),
      trueProp2Raw=matrix(rnorm(3*nrow(cpz),0,1),nrow=nrow(cpz)),
      nCells1=matrix(rnorm(3*nrow(cpz),5000,750),nrow=nrow(cpz)),
      nCells2=matrix(rnorm(3*nrow(cpz),5000,750),nrow=nrow(cpz))
    )
      #//propPlusBack1[ii,]=exp(trueProp1Raw[ii,]*sigma+to_row_vector(prop1Base[ii]))+exp(back1[ii,]);
      #//prop1Base[ii]=baseMu[ii]+repEffectPlusZero+prop1BaseRaw[ii]*baseSigma;
    warning(paste(c(exp(zz$trueProp1Raw[1,]*zz$sigma+zz$prop1BaseRaw[[1]])+exp(zz$back1[1,])),collapse=' '))
    out
  }

  selectFits2<-lapply(pairs,function(pair)compareAlleles(mod4,pair[1],pair[2],cpz))
  names(selectFits2)<-sapply(pairs,paste,collapse='-')
}

xlims<-exp(range(unlist(lapply(selectFits[names(selectFits) %in% sameScale],findLims))))
xlims2<-exp(range(unlist(lapply(selectFits[!names(selectFits) %in% sameScale],findLims))))
dummy<-lapply(names(selectFits),function(xx){message(xx);print(selectFits[[xx]],pars=c('metaFoldChange','foldChange','repEffect'))})
#pdf('cpzFits.pdf',width=4,height=4);par(mar=c(3.5,6,2,.4));lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(cpz),main=sub('-','/',xx),cols=envCols,xlims=xlims));dev.off()
pdf('out/cpzFits.pdf',width=4,height=4);par(mar=c(3.5,6,2,.4));lapply(names(selectFits),function(xx)plotFit(selectFits[[xx]],rownames(cpz),main=sub('-','/',xx),cols=NULL,xlims=if(xx %in% sameScale)xlims else xlims2));dev.off()
pdf('out/cpzRaw.pdf',width=8,height=6);par(mar=c(3,4,1,.2));lapply(pairs,function(xx)plotRaw(cpz[,grep(xx[1],colnames(cpz))],cpz[,grep(xx[2],colnames(cpz))],cols=envCols));dev.off()
xlims<-exp(range(unlist(lapply(selectFits2[names(selectFits2) %in% sameScale],findLims))))
xlims[1]<-max(xlims[1],1e-5)
xlims2<-exp(range(unlist(lapply(selectFits2[!names(selectFits2) %in% sameScale],findLims))))
pdf('out/cpzFits2.pdf',width=4,height=4);par(mar=c(3.5,6,2,.4));lapply(names(selectFits2),function(xx)plotFit(selectFits2[[xx]],rownames(cpz),main=sub('-','/',xx),cols=NULL,xlims=if(xx %in% sameScale)xlims else xlims2));dev.off()


apply(as.matrix(selectFits2[['Human-QQNVP']])[,c('metaFoldChange','metaSd')],2,mean)
apply(as.matrix(selectFits2[['QQNVP-QQNVT']])[,c('metaFoldChange','metaSd')],2,mean)
apply(as.matrix(selectFits2[['Human-QQNVT']])[,c('metaFoldChange','metaSd')],2,mean)

stats<-lapply(names(selectFits2),function(xx){
  out<-assignNames(pullRanges(selectFits2[[xx]],convertFunc=exp),rownames(cpz))
  data.frame('comparison'=xx,'stat'=rownames(out),out)
})
names(stats)<-names(selectFits)
write.csv(do.call(rbind,stats),'out/cpzStats.csv')
