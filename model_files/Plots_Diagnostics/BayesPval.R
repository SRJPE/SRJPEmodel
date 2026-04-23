rm(list = ls())
library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel) # and any others
library(rstan)
library("brms")


logit<-function(x){log(x/(1-x))}
inv_logit<-function(x){exp(x)/(1+exp(x))}

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)


par(mfcol=c(2,2),mai=c(0.75,0.75,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))

for(ii in 1:3){
  IsMain=switch(ii,F,T,T)
  MainSite=switch(ii,NA,"knights landing","tisdale")
  rname=switch(ii,"Tributary Model","Knights Landing","Tisdale")

  pCap_inputs <- prepare_pCap_inputs(mainstem=IsMain,mainstem_site=MainSite)

  Nmr=pCap_inputs$inputs$data$Nmr
  obsRec=pCap_inputs$inputs$data$Recaptures
  obsRel=pCap_inputs$inputs$data$Releases
  flow=pCap_inputs$inputs$data$mr_flow
  ind_yr=pCap_inputs$inputs$data$ind_yr


  if(IsMain==F){
    ind_trib=pCap_inputs$inputs$data$ind_trib
    load("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_trib.Rdata")
  } else {
    ind_trib=rep(1,Nmr)
    load(paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_mainstem_skew_re_",MainSite,".Rdata"))
  }
  dp=as.data.frame(pcap,pars=c("b0_pCap","b_flow","yr_re","logit_pCap","pro_sd_P"))

  print(c(rname,Nmr))


  disc.obs=rep(0,Nmr);disc.sim=disc.obs
  for(i in 1:Nmr){

    vnm=paste0("logit_pCap[",i,"]");icol=which(names(dp)==vnm);logit_pCap=dp[,icol]
    pCap=inv_logit(logit_pCap)

    e.rec=pCap*obsRel[i] #Ntrials expected recoveries

    #observed discrepancies based on Freeman and Tukey statistics
    #𝐷(𝑥;𝜃)=Σ(√𝑜𝑏𝑠−√𝑒𝑥𝑝)^2 #looks like sum the squared differences across Nmr efficiency trials. Confirmed via Gemini R code
    #The Freeman-Tukey (FT) statistic is the sum of squared, transformed values used to compute a goodness-of-fit statistic
    #disc.obs[i]=sum(sqrt(e.rec) - sqrt(obsRec[i]))^2  #sum the differences then square like Carl's code but in aaply
    disc.obs[i]=sum((sqrt(e.rec) - sqrt(obsRec[i]))^2) #sum the squared differences as in gemini example, what formula looks like, and description above

    sim.rec=rbinom(n=length(pCap),size=obsRel[i],prob=pCap) #simulated discrepencies. replace observed with simulated

    disc.sim[i]=sum((sqrt(e.rec) - sqrt(sim.rec))^2)
  }

  #Test statistic from simulated data at least as extreme as observed discrepancy
  p.value=mean(disc.sim>disc.obs) #proportion of Nmr cases when simulated discrep > observed discrep (at leastobserved discrepency > simulated discrepancy
  plot(disc.obs,disc.sim,main=paste(rname,round(p.value,digits=2),sep=" - "),bty='l')
  abline(coef=c(0,1),lty=1)
  #hist(disc.obs-disc.sim,main=paste(rname,round(p.value,digits=2),sep=" - "),bty='l')


  for(j in 1:Nmr){
      if(obsRec[j]<=2){
      sit=pCap_inputs$site_year_fit$site[ind_yr[j]]
      ryr=pCap_inputs$site_year_fit$run_year[ind_yr[j]]
      print(c(j,round(disc.obs[j],digits=0),sit,ryr,obsRel[j],obsRec[j]))
      points(disc.obs[j],disc.sim[j],pch=1,col="red")
    }
  }

}
