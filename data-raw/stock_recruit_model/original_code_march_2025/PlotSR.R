library(rstan)
library(scales)
graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

"graph_let"=function(letnum){
  usr=par("usr"); inset.x=0.03*(usr[2]-usr[1]); inset.y=0.03*(usr[4]-usr[3])
  text(usr[1]+inset.x,usr[4]-inset.y,paste(letters[letnum],")",sep=""),cex=0.75,font=2)
}

OnePlot=T
OutDir="Output/"
Nscale=0.001

UseSQcovar=T #as above but based on reduced SR dataset which only uses years when all covariates were available
#CoVarNm="si_max_flow"
for(itype in 7:7){
  DoSite=switch(itype,"ubc","ubc","lcc","lcc","mill creek","deer creek","okie dam","knights landing")
  Am=switch(itype,"rd","us","rd","us","rd","hd","hd","us")
  #CoVarNm=switch(itype,"si_max_flow","si_max_flow","rr_mean_flow","si_median_flow","si_min_flow","mi_gdd_sacramento")
  CoVarNm=switch(itype,"si_weekly_max_temp_mean","si_max_flow","rr_mean_flow","rr_median_flow",
                 "si_above_13_temp_week","si_above_13_temp_week","rr_min_flow","rr_max_flow")
  #CoVarNm="null"
  
  source("GetData.R")
  
  #Load posteriors and calculate predicted recruitment across a range of spawning stock sizes (rather than just observed ones as model does)
  Snm=paste0(OutDir,"Fit_",DoSite,"_",Am,"_",CoVarNm,".Rdata")
  load(file=Snm) 
  dp=as.data.frame(fit, pars = c("alpha","beta","gamma","sd_pro"))
  
  icol=which(names(dp)=="alpha");alpha=dp[,icol]
  icol=which(names(dp)=="beta");beta=dp[,icol]
  icol=which(names(dp)=="gamma");gamma=dp[,icol]
  icol=which(names(dp)=="sd_pro");sd_pro=dp[,icol]
  
  NxS=50;xS=seq(0,max(SP)*1.05,length.out=NxS)
  stat=matrix(nrow=NxS,ncol=3)
  for(i in 1:NxS){
    pR=xS[i]*exp(alpha + beta*xS[i])*Nscale
    stat[i,]=quantile(pR,probs=c(0.025,0.5,0.975))
    stat[i,2]=mean(pR)
  }
  
  #Compare predicted and observed recruitment as a function of spawning stock
  par(mfcol=c(1,1),mai=c(1,1,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
  #par(mfrow=c(2,2),xaxs="i",mai=c(.8,.7,.1,.5), omi=c(0.25,0.25,0.5,0.25),cex.axis=0.9)
  #par(mfcol=c(2,1),mai=c(0.85,1,0.5,1),omi=c(0.1,0.1,0.1,0.1))
  
  ymax=max(c(stat[,2],R*Nscale))*1.05
  plot(SP, R*Nscale, type='p',xlim=c(0,max(xS)),ylim=c(0,ymax),pch=19,cex=0.75,bty='l',xlab="Spawner abundance",ylab="Outmigrant abundance ('000s)")
  polygon(x = c(rev(xS), xS), y = c(rev(stat[,1]),stat[,3]), col = alpha("grey", 0.5), border = NA)
  #yr=paste0("'",substr(x=as.character(Byr), start=3,stop=4))
  yr=substr(x=as.character(Byr), start=3,stop=4)
  xoff=(max(SP))*0.025
  for(iyr in 1:Nyrs){
    text(x=SP[iyr]+xoff,y=R[iyr]*Nscale,labels=yr[iyr],cex=0.75)
  }
  lines(xS,stat[,2],lty=1)
  #mtext("Spawner Abundance ('000s)",side = 1,line = -1, outer = T,cex=1.2,font=2)
#mtext("Juvenile Abundance ('000s)",side=2, las=3,outer=T,cex=1.2,font=2,line=-1)

  #Plot covariate effect on SR graph
  if(CoVarNm!="mu" & CoVarNm!="null" ){
    for(iyr in 1:Nyrs){
      pR1=mean(SP[iyr]*exp(alpha + beta*SP[iyr]))*Nscale
      pR2=mean(SP[iyr]*exp(alpha + beta*SP[iyr]+ gamma*X[iyr]))*Nscale
      if((pR2>pR1 & (R[iyr]*Nscale)>pR1) | (pR2<pR1 & (R[iyr]*Nscale)<pR1)){
        xcol="blue"
      } else {
        xcol="red"
      }
      lines(x=c(SP[iyr],SP[iyr]),y=c(pR1,pR2),lty=1,lwd=2,col=xcol)
    }
    legend("bottomright",cex=0.9,legend=c("@ avg cov value", "covariate effect (consistent)","covariate effect (inconsistent)"),lty=c(1,1),lwd=c(2,2),col=c("black","blue","red"),bty='n')
  }
  #graph_let(1)
  #mtext("Spawner Abundance",side=1, line=-1,las=1,outer=T,cex=1.3,font=2)
  #mtext("Outmigrant abundance ('000s)",side=2, line=-1,las=3,outer=T,cex=1.3,font=2)
  mtext(paste0(DoSite,"_",Am," (",CoVarNm,")"),side=3, line=-1,outer=T,cex=1.1,font=2)
    
  if(OnePlot==F){
    #Plot covariate relationshp
    if(CoVarNm!="mu" & CoVarNm!="null"){
      
      Nbins=50
      x=seq(min(X0)*1,max(X0)*1,length.out=Nbins)
      xstd=(x-mean(X0))/sd(X0)
      
      stat=matrix(nrow=Nbins,ncol=3)
      muS=mean(SP)
      for(i in 1:Nbins){
          pred=muS*exp(alpha+beta*muS + gamma*xstd[i])*Nscale
          stat[i,]=quantile(pred,probs=c(0.025,0.5,0.975))
          stat[i,2]=mean(pred)
      }
      ymax=max(c(stat[,2],R*Nscale))*1.25
      plot(X0,R*Nscale,pch=19,cex=0.75,xlim=c(min(x),max(x)),ylim=c(0,ymax),bty='l',xlab=CoVarNm,ylab="")
      polygon(x = c(rev(x), x), y = c(rev(stat[,1]), stat[,3]), col = alpha("grey", 0.25), border = NA)
      lines(x,stat[,2],lty=1,lwd=2)
      xoff=(max(X0)-min(X0))*0.025
      for(iyr in 1:Nyrs){
        text(x=X0[iyr]+xoff,y=R[iyr]*Nscale,labels=yr[iyr],cex=0.75)
        p1=mean(muS*exp(alpha+beta*muS + gamma*X[iyr]))*Nscale
        p2=mean(SP[iyr]*exp(alpha+beta*SP[iyr] + gamma*X[iyr]))*Nscale
        lines(x=c(X0[iyr],X0[iyr]),y=c(p1,p2),col="red")
      }
      legend("topright",legend=c("@ avg spawners", "with spawner effect"),lty=c(1,1),lwd=c(2,2),col=c("black","red"),bty='n')
      graph_let(2)
      mtext("Outmigrant abundance ('000s)",side=2, line=-1,las=3,outer=T,cex=1.3,font=2)
      mtext(paste0(DoSite,"_",Am),side=3, line=-1,outer=T,cex=1.3,font=2)
      
      #Stock-recruit at lowest and highest covariate values
      stat=array(dim=c(2,NxS,3))
      for(i in 1:NxS){
        pred=xS[i]*exp(alpha + beta*xS[i] + gamma*min(X))*Nscale
        stat[1,i,]=quantile(pred,probs=c(0.025,0.5,0.975));  stat[1,i,2]=mean(pred)
        pred=xS[i]*exp(alpha + beta*xS[i] + gamma*max(X))*Nscale
        stat[2,i,]=quantile(pred,probs=c(0.025,0.5,0.975));stat[2,i,2]=mean(pred)
      }
      plot(xS,stat[1,,2],type='l',lwd=2,col="blue",bty='n',ylim=c(0,max(stat[,,3])),xlab="Spawner Abundance",ylab="")
      polygon(x = c(rev(xS), xS), y = c(rev(stat[1,,1]), stat[1,,3]), col = alpha("blue", 0.25), border = NA)
      lines(xS,stat[2,,2],col="red")
      polygon(x = c(rev(xS), xS), y = c(rev(stat[2,,1]), stat[2,,3]), col = alpha("red", 0.25), border = NA)
      legend("topleft",legend=c("min cov value", "max cov value"),lty=c(1,1),lwd=c(2,2),col=c("blue","red"),bty='n')
      graph_let(3)
    }
    
    #par(mfcol=c(2,1),mai=c(0.85,1,0.25,1),omi=c(0.1,0.1,0.1,0.1))
    
    #Juveniles with error time series and covariates on second axis
    statR=matrix(nrow=Nyrs,ncol=2)
    for(iyr in 1:Nyrs){
      simR=exp(rnorm(n=1000,mean=mu_obslgR[iyr],sd=sd_obslgR[iyr]))*Nscale
      statR[iyr,]=quantile(simR,probs=c(0.025,0.975))
    }
    
    x=barplot(R*Nscale,space=0.1,names.arg="",col="grey",bty='l',xlab="Brood Year",ylab="",axes=F,ylim=c(0,max(statR)),xlim=c(0,Nyrs+1))
    axis(2);axis(1,at=x,labels=yr,las=3,cex.axis=0.9)
    abline(h=mean(R)*0.001,lty=2)
    for(iyr in 1:Nyrs){
      arrows(x0=x[iyr],x1=x[iyr],y0=statR[iyr,1],y1=statR[iyr,2],angle=90,code=3,length=0.025)
    }
    
    if(CoVarNm!="mu" & CoVarNm!="null"){
      par(new=T)
      plot(x,X0,pch=19,cex=1.3,axes=F,ylab="",xlab="",ylim=c(min(X0)*0.95,max(X0)*1.15),xlim=c(0,Nyrs+1)); axis(4)
      mtext(CoVarNm, side=4, line=3,las=3)
      legend("top",legend=c("outmigrants","covariate"),pch=c(15,19),col=c("grey","black"),bty='y',ncol=2)
    }
    graph_let(4)
    
    mtext("Outmigrant abundance ('000s)",side=2, line=-1,las=3,outer=F,cex=1.3,font=2)
    mtext(paste0(DoSite,"_",Am),side=3, line=-1,outer=F,cex=1.3,font=2)
  }
  #Do a forecast assuming spawning stock is hitorical average spawner abundance, and covariate is at average historical conditon
  #Differentiate between parameter uncertainty and process error in plot
  par(mfcol=c(1,1),mai=c(1,1,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
  muSP=mean(SP)
  pred1=muSP*exp(alpha + beta*muSP)*Nscale
  nsims=length(pred1)
  dev=rnorm(n=nsims,mean=0,sd=mean(sd_pro))
  pred2=pred1*exp(dev)
  print("")
  print(c(DoSite," ",Am))
  print(quantile(pred1,prob=c(0.025,0.5,0.975)));print(sd(pred1)/mean(pred1))
  print(quantile(pred2,prob=c(0.025,0.5,0.975)));print(sd(pred2)/mean(pred2))
}