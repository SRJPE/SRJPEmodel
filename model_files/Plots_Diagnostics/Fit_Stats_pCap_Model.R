#Get fit statistics for pCap model. This includes r2, r2 with predictions based on fixed effects only, and bayesian p-value

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

  if(IsMain==F){
    pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
  } else {
    pCap_inputs <- prepare_pCap_inputs(model_type = "one_site", skew = T, site_selection = MainSite)
  }

  Nmr=pCap_inputs$inputs$data$Nmr
  obsRec=pCap_inputs$inputs$data$Recaptures
  obsRel=pCap_inputs$inputs$data$Releases
  flow=pCap_inputs$inputs$data$mr_flow
  ind_yr=pCap_inputs$inputs$data$ind_yr


  if(IsMain==F){
    ind_trib=pCap_inputs$inputs$data$ind_trib
    load("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_all_sites.Rdata")
  } else {
    ind_trib=rep(1,Nmr)
    load(paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_one_site_skew_re_",MainSite,".Rdata"))
  }
  dp=as.data.frame(pcap,pars=c("b0_pCap","b_flow","yr_re","logit_pCap","pro_sd_P"))
  Nsims=dim(dp)[1]

  print(c(rname,Nmr))

  pred_pCap=vector(length=Nmr);obs_pCap=pred_pCap;pred_pCap_fixed=pred_pCap;pred_pCap_fixed2=pred_pCap
  disc.obs=rep(0,Nmr);disc.sim=disc.obs
  est_fe_re=matrix(nrow=Nsims,ncol=Nmr);est_re=est_fe_re

  for(i in 1:Nmr){

    vnm=paste0("logit_pCap[",i,"]");icol=which(names(dp)==vnm);logit_pCap=dp[,icol]

    pCap=inv_logit(logit_pCap)

    e.rec=pCap*obsRel[i] #Ntrials expected recoveries

    #observed discrepancies based on Freeman and Tukey statistics
    #𝐷(𝑥;𝜃)=Σ(√𝑜𝑏𝑠−√𝑒𝑥𝑝)^2 #looks like sum the squared differences across Nmr efficiency trials. Confirmed via Gemini R code
    #The Freeman-Tukey (FT) statistic is the sum of squared, transformed values used to compute a goodness-of-fit statistic
    #disc.obs[i]=sum(sqrt(e.rec) - sqrt(obsRec[i]))^2  #sum the differences then square like Carl's code but in aaply so not obvious that sum was across mcmc trials
    #for each efficiency trial i.
    disc.obs[i]=sum((sqrt(e.rec) - sqrt(obsRec[i]))^2) #sum the squared differences as in gemini example, what formula looks like, and description above

    sim.rec=rbinom(n=Nsims,size=obsRel[i],prob=pCap) #simulated discrepencies. replace observed with simulated

    disc.sim[i]=sum((sqrt(e.rec) - sqrt(sim.rec))^2)

    pred_pCap[i]=mean(pCap) #with prodcess error
    obs_pCap[i]=obsRec[i]/obsRel[i]

    #Compute pred_pCap based on fixed and year effects only
    if(IsMain==T){
      vnm="b0_pCap";icol1=which(names(dp)==vnm)
      vnm="b_flow";icol2=which(names(dp)==vnm)
      vnm=paste0("yr_re[",ind_yr[i],"]");icol3=which(names(dp)==vnm)
    } else {
      vnm=paste0("b0_pCap[",ind_trib[i],"]");icol1=which(names(dp)==vnm)
      vnm=paste0("b_flow[",ind_trib[i],"]");icol2=which(names(dp)==vnm)
      vnm=paste0("yr_re[",ind_yr[i],"]");icol3=which(names(dp)==vnm)
    }
    pred_pCap_fixed[i]=mean(inv_logit(dp[,icol1] + dp[,icol2]*flow[i] + dp[,icol3]))
    pred_pCap_fixed2[i]=mean(inv_logit(dp[,icol1] + dp[,icol2]*flow[i]))

    #For Gelman and Pardoe (2006) pR2 statistic
    est_fe_re[,i]=logit_pCap #fixed + random effects
    #backout prod_dev_P since it wasn't saved.  logit_pCap=fixed effects + pro_dev_P, so pro_dev_p = logit_pCap-fixed effects
    est_re[,i]=logit_pCap-(dp[,icol1] + dp[,icol2]*flow[i] + dp[,icol3]) #random effects only
  }

  #Test statistic from simulated data at least as extreme as observed discrepancy
  p.value=round(mean(disc.sim>disc.obs),digits=2) #proportion of Nmr cases when simulated discrep > observed discrep (at leastobserved discrepency > simulated discrepancy
  plot(disc.obs,disc.sim,main=paste(rname, p.value,sep=" - "),bty='l')
  abline(coef=c(0,1),lty=1)
  #hist(disc.obs-disc.sim,main=paste(rname,round(p.value,digits=2),sep=" - "),bty='l')

  for(j in 1:Nmr){
      if(obsRec[j]<=2){
      sit=pCap_inputs$site_year_fit$site[ind_yr[j]]
      ryr=pCap_inputs$site_year_fit$run_year[ind_yr[j]]
      #print(c(j,round(disc.obs[j],digits=0),sit,ryr,obsRel[j],obsRec[j]))
      points(disc.obs[j],disc.sim[j],pch=1,col="red")
    }
  }


  #Complte multi-level r2 calcs
  Var_fe_re=vector(length=Nsims);Var_re=Var_fe_re
  mure=colMeans(est_re) #the mean of each random effect across posterior samples
  Var_mure=var(mure) #the variance across the posterior means of each random effect

  for(isim in 1:Nsims){
    Var_fe_re[isim]=var(est_fe_re[isim,]) #the variance across logit-transformed survival, movement, or pCap for each posterior sample
    Var_re[isim]=var(est_re[isim,])
  }
  pR2=round(1-mean(Var_re)/mean(Var_fe_re),digits=3) #eqn. 6 (proportion of variation in rate explained by fixed effects)
  lambda=round(1-Var_mure/mean(Var_re),digits=3)     #eqn 7 (extent of pooling of residuals)

  r2=round(cor(pred_pCap,obs_pCap)^2,digits=3)
  r2_fixed=round(cor(pred_pCap_fixed,obs_pCap)^2,digits=3)
  r2_fixed2=round(cor(pred_pCap_fixed2,obs_pCap)^2,digits=3)

  print(c(rname, p.value, r2, r2_fixed ,r2_fixed2, pR2, lambda))

}
