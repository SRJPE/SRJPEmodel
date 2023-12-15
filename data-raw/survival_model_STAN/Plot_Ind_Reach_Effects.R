library("rstan")
library("scales")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)


source("GetData.R")

par(mfcol=c(2,2),mai=c(.7,0.6,0.2,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1,font.main=1,cex.lab=1,font.lab=2)

varnm="SurvWoodSac"
for(ivar in 1:3){
  ext=switch(ivar,"Q","Vel","Wt")
  load(file=paste0("Results/fit_CovInd_Reach_",ext,".Rdata"))
  
  survmat=as.data.frame(fit, pars = c(varnm))
  stats=matrix(nrow=Nseq,ncol=3)
  for(ix in 1:Nseq){
      vname=paste0(varnm,"[",ix,"]")
      icol=which(names(survmat)==vname)
      stats[ix,]=quantile(survmat[,icol],prob=c(0.025,0.5,0.975))
      stats[ix,2]=mean(survmat[,icol])
  }
  
  if(ivar==1){
    X=cbind(d$Flow_Releasepoint,d$Flow_Woodson.Br, d$Flow_ButteBridge,d$Flow_Sacramento)  #Flow
    xlabel="Discharge (m3/sec)"
  } else if (ivar==2){
    X=cbind(d$Vel_Releasepoint,d$Vel_Woodson.Br, d$Vel_ButteBridge,d$Vel_Sacramento) #Velocity
    xlabel="Velocity (m/sec)"
  } else if (ivar==3){
    X=cbind(d$Temp_Releasepoint,d$Temp_Woodson.Br, d$Temp_ButteBridge,d$Temp_Sacramento) #Temperature
    xlabel="Water Temperature (C)"
  }
  Xstd=matrix(nrow=Nind,ncol=Nreaches); muX=mean(X);sdX=sd(X)
  for(j in 1:Nreaches){Xstd[,j]=(X[,j]-muX)/sdX} #Calculate the standardized values used to predict covariate function values in model
  Xseq_std=seq(from=min(Xstd),to=max(Xstd),length.out=Nseq) 
  Xseq=Xseq_std*sdX+muX #transform back to original units
  
  #RE_sd=round(mean(unlist(as.data.frame(fit, pars = c("RE_sd")))),digits=2)
  #RE_sd_sd=round(sd(unlist(as.data.frame(fit, pars = c("RE_sd")))),digits=2)
  plot(Xseq,stats[,2],type='l',lwd=2,bty='n',xlab=xlabel,ylab="",ylim=c(0,0.6)) #max(stats),main=paste0("RE_sd = ",RE_sd," (", RE_sd_sd,")")
  polygon(x = c(rev(Xseq),Xseq), y = c(rev(stats[,1]), stats[,3]), col = alpha("grey", 0.5), border = NA)
  abline(h=seq(0.1,0.6-0.1,by=0.1),col="grey",lty=2)
}
mtext("Woodson-Sacramento Survival Rate",side=2, line=1,las=3,outer=T,cex=1.2,font=2)
    