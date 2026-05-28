rm(list=ls(all=TRUE))
library(rstan)
library(scales)
library(dplyr)
library(brms)
library(SRJPEdata)
library(SRJPEmodel)

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,height=12,width=16);par(las=1)

outdir="C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/0p5minP"
#pdf(file=paste0(outdir,"/SiteYear_0p5minP.pdf"))

lb=0.025;ub=0.975	#0.975
Nscale=0.001		#abundance is in units of '000s of fish

d0a=SRJPEdata::years_to_include_rst_data

d0=subset(d0a,site!="lbc" & site!="red bluff diversion dam")

siteO=c("ucc","lcc","ubc","mill creek","deer creek","okie dam","steep riffle","gateway riffle","eye riffle","herringer riffle",
        "live oak","sunset pumps","hallwood","tisdale","knights landing")

uniqsite=siteO
Nsites=length(uniqsite)

for(isite in 1:Nsites){#Nsites

  d1=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]

  if(DoSite=="ucc"){#load pcap fit object and pcap inputs
    load(file="C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_all_sites.rdata")
    model_type="all_sites"
    if(isite==1) pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
  } else if (DoSite=="tisdale" | DoSite=="knights landing"){
    load(file=paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_one_site_skew_re_",DoSite,".rdata"))
    model_type="one_site_skew"#_re
    pCap_inputs <- prepare_pCap_inputs(model_type = "one_site_skew", skew = T, site_selection = DoSite)

  }

  uniqyr=unique(sort(d1$run_year))
  Nyrs=length(uniqyr)

  for(iyr in 1:Nyrs){#1:Nyrs
  #for(iyr in Nyrs:Nyrs){#1:Nyrs

    DoYr=uniqyr[iyr]
    RunNm=paste(DoSite,DoYr,sep="-")

    abundance_inputs <- prepare_abundance_inputs(site = DoSite,
                                                 run_year = DoYr,
                                                 pCap_model_type = model_type,
                                                 min_pCap_mult = 0.5,
                                                 pCap_model_object = pcap)


    print(c(DoSite,DoYr,abundance_inputs$run_year,abundance_inputs$run_year_id))

    Nstrata=abundance_inputs$inputs$data$Nstrata
    Trib=abundance_inputs$sites_fit
    Ntribs=length(Trib)
    Jdt=abundance_inputs$week_date

    fn_post=paste0(outdir,"/",DoSite,"_",DoYr,".rdata")
    load(file=fn_post)

    par(mfcol=c(2,1),mai=c(0.8,1,0.2,0.7), omi=c(0.1,0.1,0.3,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)

    ####Plot weekly abundance ####################
    dp=abundance$sims.list$Ntot
    Ntot=as.double(quantile(dp,prob=c(lb,0.5,ub)))*Nscale
    Ntot_cv=round(100*sd(dp)/mean(dp),digits=0)
    NtotLab=paste0(round(Ntot[2],digits=0)," (",round(Ntot[1],digits=0)," - ", round(Ntot[3],digits=0), ") cv=",Ntot_cv,"%")

    #Vnm="N";dp=as.data.frame(abundance,pars=Vnm)
    dp=as.data.frame(abundance$sims.list$N)
    N=matrix(nrow=Nstrata,ncol=3)
    for(i in 1:Nstrata){
      N[i,]=as.double(quantile(dp[,i]*Nscale,prob=c(lb,0.5,ub)));N[i,2]=mean(dp[,i])*Nscale
    }


    u=abundance_inputs$inputs$data$u
    Uwc_ind=abundance_inputs$inputs$data$Uwc_ind
    #irecs_N=abundance_inputs$lt_pCap_U_data$Uind_wMR
    irecs_N=which(!is.na(abundance_inputs$lp_data$number_released) &
                    !is.na(abundance_inputs$lp_data$number_recaptured))

    pN=0;Releases=NA;Recaptures=NA
    if(length(irecs_N)>0){
      Releases=abundance_inputs$lp_data$number_released[irecs_N]
      Recaptures=abundance_inputs$lp_data$number_recaptured[irecs_N]

      irecs_u=which(is.na(match(Uwc_ind,irecs_N))==F)
      #pN=(u[irecs_u]/((Recaptures)/Releases))*Nscale
      pN=(u[irecs_u]/((Recaptures+0.01)/Releases))*Nscale
      #sdN=sqrt(Nscale^2*(as.numeric(Releases)*as.numeric(Recaptures)*(as.numeric(Releases)-as.numeric(Recaptures))*u[irecs_N])/(Recaptures^2*Recaptures))
      ibad=which(Recaptures==0);  pN[ibad]=NA
    }

    #Jwk=abundance_inputs$weeks_fit
    xmain=paste0(RunNm," Ntot=",NtotLab)
    ylims=c(-0.1,max(c(max(N),max(pN,na.rm=T)))*1.1)
    bcols=rep("grey",Nstrata)
    for(i in 1:Nstrata){if(is.na(match(i,Uwc_ind))==T) bcols[i]="red"} #bars with no catch = no sampling =red.

    x=barplot(height=N[,2],col=bcols,names.arg=Jdt,cex.names=0.7, las=3,axis.lty=1,ylim=ylims,  bty='l',main=xmain, xlab="",ylab="Abundance ('000s)")
    #points(x=x[round(kknots,digits=0)],y=rep(ylims[2]*0.8,Nknots),pch=25,bg="green",cex=1.5)#knot position
    k=0
    for(i in 1:Nstrata){#error bars and catch
      arrows(x0=x[i],x1=x[i],y0=N[i,1],y1=N[i,3],angle=90,code=3,length=0.025,col="black")
      if(is.na(match(i,Uwc_ind))==T){
        points(x[i],ylims[2]*0.9,pch=19,col="red",cex=0.8)#identify strata with no sampling (red bars may not show up if estimated abundance is very low)
      } else {
        k=k+1
        par(srt=90);text(x=x[i],y=ylims[2]*0.9,labels=u[k],cex=0.7)
        text(x=x[i],y=ylims[2]*0.7,labels=round(abundance_inputs$inputs$data$effort[i],digits=2),cex=0.7)
        par(srt=0)
      }
    }

    #Plot Peterson estimates and CL's
    coff=0.1
    nn=length(irecs_N)
    if(nn>0){
      for(i in 1:nn){
        points(x[irecs_N[i]]-coff,pN[i],pch=21,cex=1.5,col='blue')	#data-driven predictions of N
      }
    }
    legend("right",legend=c("Bayesian","Peterson","No sampling or no releases","Discharge"),pch=c(22,21,19,NA),pt.cex=c(rep(1.2,3),NA),pt.bg=c("grey",NA,NA,NA),col=c("grey","blue","red","darkgreen"),lwd=c(NA,NA,NA,1.2),cex=0.8,bty='n')

    Flow=abundance_inputs$catch_flow/1000;xlab="Discharge (kcfs)"
    #Flow=inputs$catch_flow_raw/1000;xlab="Discharge (kcfs)"
    #df0=weekly_juvenile_abundance_catch_data
    #df1=subset(df0,site==DoSite & run_year==DoYr) #would still be out of order and needs to do for one lifestage only
    par(new=T)
    plot(x-coff,Flow,ylim=c(min(Flow,na.rm=T),max(Flow,na.rm=T)*1.2), xlim=c(0,trunc(max(x))+1),type='o',pch=19,cex=0.7,lwd=2,col="darkgreen",axes=F,xlab="",ylab="");axis(4,cex.axis=0.9)
    mtext(xlab,side=4, las=3,adj=0.5,outer=T,cex=1.0,font=2,line=-1,col="darkgreen")



    ###Plot weekly capture probability #########################
    dp=abundance$sims.list$pCap_U
    pCap=matrix(nrow=Nstrata,ncol=3)
    for(i in 1:Nstrata){

      #identify records where pCap estimate is equal or above assumed minimum used to set the upper bound on N (lgN_max[iwk])
      #Only select these posterior samples for the plot, since that is what is effectively used in the model based on upper contraint on N
      irecs=which(dp[,i]>=abundance_inputs$min_pCap)
      if(length(irecs)==0) irecs=1:dim(dp)[1]
      pCap[i,]=as.double(quantile(dp[irecs,i],prob=c(lb,0.5,ub))); pCap[i,2]=mean(dp[irecs,i])
      #pCap[i,]=as.double(quantile(dp[,i],prob=c(lb,0.5,ub))); pCap[i,2]=mean(dp[,i])
    }
    #ylims=c(0,max(pCap,na.rm=T)*1.1)
    ylims=c(0,max(pCap)*1.1)

    bcols=rep("red",Nstrata);bcols[irecs_N]="grey"
    x=barplot(height=pCap[,2],col=bcols,names.arg=Jdt,cex.names=0.7, las=3,axis.lty=1,ylim=ylims,  bty='l', xlab="First Date of Week",ylab="Capture Probability")

    irecs=which(is.na(Releases)==F)
    if(length(irecs)>0) points(x[irecs_N]-coff,Recaptures/Releases,pch=21,cex=1.5,col='blue')

    k=0
    for(i in 1:Nstrata) {
      arrows(x0=x[i],x1=x[i],y0=pCap[i,1],y1=pCap[i,3],angle=90,code=3,length=0.025,col="black")
      if(length(which(irecs_N==i))>0){
        par(srt=90)
        k=k+1
        text(x=x[i],y=ylims[2]*0.95,labels=Recaptures[k],cex=0.7)
        text(x=x[i],y=ylims[2]*0.82,labels=Releases[k],cex=0.7)
        par(srt=0)
      }
    }
    par(new=T)
    plot(x-coff,Flow,ylim=c(min(Flow,na.rm=T),max(Flow,na.rm=T)*1.2), xlim=c(0,trunc(max(x))+1),type='o',pch=19,cex=0.7,lwd=2,col="darkgreen",axes=F,xlab="",ylab="");axis(4,cex.axis=0.9)

  }#iyr
}#isite
#dev.off()#if writing plots to a file.



