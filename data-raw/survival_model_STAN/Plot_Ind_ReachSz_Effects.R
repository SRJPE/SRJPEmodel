library("rstan")
library("scales")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)


source("GetData.R")

par(mfcol=c(2,2),mai=c(.7,0.6,0.2,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1,font.main=1,cex.lab=1,font.lab=2)

varnm="SurvWoodSacSz"
for(ivar in 1:3){
  ext=switch(ivar,"Vel_FL","Vel_WGT","Vel_CF")
  load(file=paste0("Results/fit_CovInd_ReachSz_",ext,".Rdata"))
  
  survmat=as.data.frame(fit, pars = c(varnm))
  stats=matrix(nrow=Nseq,ncol=3)
  for(ix in 1:Nseq){
    vname=paste0(varnm,"[",ix,"]")
    icol=which(names(survmat)==vname)
    stats[ix,]=quantile(survmat[,icol],prob=c(0.025,0.5,0.975))
    stats[ix,2]=mean(survmat[,icol])
  }
  
  if(ivar==1){
    xlabel="Fork Length (mm)"
    Xsz=FL
    xlims=c(70,140)
  } else if (ivar==2){
    Xsz=WGT
    xlims=c(5,30)
    xlabel="Weight (g)"
  } else if (ivar==3){
    Xsz=CF*1000
    xlims=c(0.8,1.8)
    xlabel="Condition Factor ('000s)"
  }
  muXsz=mean(Xsz);sdXsz=sd(Xsz)
  Xsz_std=(Xsz-muXsz)/sdXsz
  Xsz_std_seq=seq(from=min(Xsz_std),to=max(Xsz_std),length.out=Nseq)
  Xsz_seq=Xsz_std_seq*sdXsz+muXsz #transform back to original units
  
  plot(Xsz_seq,stats[,2],type='l',lwd=2,bty='n',xlab=xlabel,ylab="",xlim=xlims,ylim=c(0,0.6),axes=F) #max(stats),main=paste0("RE_sd = ",RE_sd," (", RE_sd_sd,")")
  axis(2);axis(1,las=3)
  polygon(x = c(rev(Xsz_seq),Xsz_seq), y = c(rev(stats[,1]), stats[,3]), col = alpha("grey", 0.5), border = NA)
  abline(h=seq(0.1,0.6-0.1,by=0.1),col="grey",lty=2)
}
mtext("Woodson-Sacramento Survival Rate",side=2, line=1,las=3,outer=T,cex=1.2,font=2)
