#Compare annual abundance estimates based on different lgN_max priors

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,height=12,width=16);par(las=1)


#lb=0.025;ub=0.975	#0.975
lb=0.1;ub=0.9
Nscale=0.001		#abundance is in units of '000s of fish

par(mfrow=c(4,4),mai=c(0.5,0.6,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)


ModNm=c("0p1minP","0p5minP")
outdir=c(paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/",ModNm[1]),
         paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/",ModNm[2]))

d0=subset(d0a,site!="lbc" & site!="red bluff diversion dam")

siteO=c("ucc","lcc","ubc","mill creek","deer creek","okie dam","steep riffle","gateway riffle","eye riffle","herringer riffle",
        "live oak","sunset pumps","hallwood","tisdale","knights landing")

uniqsite=siteO
Nsites=length(uniqsite)

for(isite in 1:Nsites){

  d1=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]
  uniqyr=unique(sort(d1$run_year))
  Nyrs=length(uniqyr)

  #get all files with where DoSite is found in name
  file.ls=cbind(list.files(path=outdir[1],pattern=DoSite),list.files(path=outdir[1],pattern=DoSite))

  Ntot=array(dim=c(2,Nyrs,3));Ncv=matrix(nrow=2,ncol=Nyrs)
  for(iyr in 1:Nyrs){

    DoYr=uniqyr[iyr]
    RunNm=paste(DoSite,DoYr,sep="-")

    for(imod in 1:2){
      load(file=paste0(outdir[imod],"/",DoSite,"_",DoYr,".rdata"))
      dp=abundance$sims.list$Ntot
      Ntot[imod,iyr,]=as.double(quantile(dp,prob=c(lb,0.5,ub)))*Nscale
      Ncv[imod,iyr]=sd(dp)/mean(dp)
    }
  }#iyr

  if(iyr==Nyrs){ #to inspect min_pCap used to calculate lg_Nmax prior or pCap for site
    if(DoSite=="ucc"){ #load pcap fit object needed by prepare_abundance_inputs. Sites ordered by trib then mainstem so only load trib pcap object for first trib
      IsMain=F
      pCap_inputs <- prepare_pCap_inputs(mainstem=IsMain)
      #load(file="C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_trib.rdata")

    } else if (DoSite=="tisdale" | DoSite=="knights landing"){
      IsMain=T
      pCap_inputs <- prepare_pCap_inputs(mainstem=IsMain,mainstem_site=DoSite)
      irecs=1:length(pCap_inputs$inputs$data$Recaptures)
      #load(file=paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_mainstem_skew_re_",DoSite,".rdata"))
    }
    if(IsMain==F) {
      itrib=which(pCap_inputs$sites_fit==DoSite)
      irecs=which(pCap_inputs$inputs$data$ind_trib==itrib)

    }
    mu_pCap=mean(pCap_inputs$inputs$data$Recaptures[irecs]/pCap_inputs$inputs$data$Releases[irecs])
    #abundance_inputs <- prepare_abundance_inputs(site = DoSite,run_year = DoYr, effort_adjust=T, pcap_model_object = pcap)
    #min_pCap=abundance_inputs$min_pCap #the lowest observed non-zero pCap * multiplier (e.g., 1, 0.5, 0.1): what is used to define lgN_max
  }

  ylims=c(0,max(Ntot))
  xmain=paste0(DoSite," (",round(mu_pCap*100,digits=2),"%  ",round(mean(Ntot[1,,2])/mean(Ntot[2,,2]),digits=2),"  ", round(mean(Ncv[1,])/mean(Ncv[2,]),digits=2),")")
  if(isite==Nsites & iyr==Nyrs){
    x=barplot(height=Ntot[,,2],beside=T,names.arg=uniqyr,cex.names=0.7, las=3,axis.lty=1,ylim=ylims, bty='l',main=xmain,legend.text=ModNm,args.legend=c(bty='n'))
  } else {
    x=barplot(height=Ntot[,,2],beside=T,names.arg=uniqyr,cex.names=0.7, las=3,axis.lty=1,ylim=ylims, bty='l',main=xmain)
  }
  for(iyr in 1:Nyrs){
    for(imod in 1:2){
      arrows(x0=x[imod,iyr],x1=x[imod,iyr],y0=Ntot[imod,iyr,1],y1=Ntot[imod,iyr,3],angle=90,code=3,length=0.025,lwd=1)
    }
  }

}#isite
mtext("Run Year",side = 1,line = -1, outer = T,cex=1.1,font=2)
mtext("Abundance ('000s)",side=2, las=3,outer=T,cex=1.1,font=2,line=-1)
