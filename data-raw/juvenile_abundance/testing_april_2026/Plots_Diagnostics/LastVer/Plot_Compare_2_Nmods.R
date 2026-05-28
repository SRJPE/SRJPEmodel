#one-one plot of annual abundance estimates from two altnerate priors on lgN_max for abundance model


graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,height=12,width=16);par(las=1)



lb=0.025;ub=0.975	#0.975
Nscale=0.001		#abundance is in units of '000s of fish

RunMainstem=F

if(RunMainstem==F){
  d0=subset(d0a,site!="knights landing"& site!="tisdale" & site!="red bluff diversion dam")

  outdir=c("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/trib",
           "C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/trib_minP")

  par(mfrow=c(4,4),mai=c(0.5,0.6,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
} else {
  d0=subset(d0a,site=="knights landing" | site=="tisdale")
  outdir=c("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/mainstem",
           "C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/mainstem_minP")

  par(mfrow=c(2,1),mai=c(1,1,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
}


uniqsite=unique(d0$site)
Nsites=length(uniqsite)

for(isite in 1:Nsites){

  d1=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]
  uniqyr=unique(sort(d1$run_year))
  Nyrs=length(uniqyr)

  #get all files with where DoSite is found in name

  file.ls=cbind(list.files(path=outdir[1],pattern=DoSite),list.files(path=outdir[1],pattern=DoSite))

  Ntot=array(dim=c(2,Nyrs,3))
  for(iyr in 1:Nyrs){#1:Nyrs

    DoYr=uniqyr[iyr]
    RunNm=paste(DoSite,DoYr,sep="-")

    for(imod in 1:2){
      load(file=paste0(outdir[imod],"/",DoSite,"_",DoYr,".rdata"))
      dp=abundance$sims.list$Ntot
      Ntot[imod,iyr,]=as.double(quantile(dp,prob=c(lb,0.5,ub)))*Nscale

    }
  }#iyr

  ylims=c(0,max(Ntot))
  x=barplot(height=Ntot[,,2],beside=T,names.arg=uniqyr,cex.names=0.7, las=3,axis.lty=1,ylim=ylims, bty='l',main=DoSite)
  for(iyr in 1:Nyrs){
    for(imod in 1:2){
      arrows(x0=x[imod,iyr],x1=x[imod,iyr],y0=Ntot[imod,iyr,1],y1=Ntot[imod,iyr,3],angle=90,code=3,length=0.025,lwd=1)
    }
  }

}#isite
mtext("Run Year",side = 1,line = -1, outer = T,cex=1.1,font=2)
mtext("Abundance ('000s)",side=2, las=3,outer=T,cex=1.1,font=2,line=-1)







