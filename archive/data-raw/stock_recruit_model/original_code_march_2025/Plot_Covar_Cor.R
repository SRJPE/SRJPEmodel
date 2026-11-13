library("corrplot")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

"graph_let"=function(letnum){
  usr=par("usr"); inset.x=0.03*(usr[2]-usr[1]); inset.y=0.03*(usr[4]-usr[3])
  text(usr[1]+inset.x,usr[4]-inset.y,paste(letters[letnum],")",sep=""),cex=0.75,font=2)
}

#Am="us"
Am="rd"

for(isite in 1:2){
  DoSite=switch(isite,"ubc","lcc")
  fn_cov=paste0(OutDir,"CovSQ_",DoSite,"_",Am,".dat")
  dCov=read.table(file=fn_cov,header=T)
  
  par(mfcol=c(1,1),mai=c(0.9,0.9,0.1,0.1),omi=c(0.1,0.1,0.1,0.1))
  
  #covariates used to predict esapement-outmigration SR relationship
  CoVarNm=unique(dCov$Variable)
  Ncovars=length(CoVarNm)
  
  dCov1=subset(dCov,Variable==CoVarNm[1])
  Nyrs=dim(dCov1)[1]

  codat=matrix(nrow=Nyrs,ncol=Ncovars)
  for(ivar in 1:Ncovars){
    dCov1=subset(dCov,Variable==CoVarNm[ivar])
    codat[,ivar]=dCov1$Value#in case there are dups as there was for lcc
  }
  codat=as.data.frame(codat)
  colnames(codat)=CoVarNm
  
  cormat=cor(codat,codat)
  corrplot(cormat,method="number",type="upper",number.cex=0.7,tl.cex=0.8)
  graph_let(isite)
  mtext(paste0(DoSite,"_",Am),side=3, line=-1,outer=T,cex=1.3,font=2)

}
