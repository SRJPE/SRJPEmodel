rm(list = ls())
library("rstan")
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
library(SRJPEdata)
library(SRJPEmodel)

nchains=3;niter=2000
RunMainstem=T

#shape and rate pars for gamma prior on sd_pro
pr_sd_pro=c(20,10)#median of 2, 95 %CI's 1.2-3 which covers range for sites where model converges without a gamma prior


for(irun in 1:4){ #1:4without and with covariate effect on phi for without and with lag-1 autocorr
  
  if(irun==1){
    ModNm="BetaDevHBMRT"
    OutDir="Output/NoCovEffect/"
    UseCovX=c(0,0)
  
  } else if (irun==2){
    ModNm="BetaDevHBMRT"
    OutDir="Output/WithCovEffect/"
    UseCovX=c(1,0)
  
  } else if(irun==3){
    ModNm="BetaDevHBMRT_lag1"    
    OutDir="Output/lg1NoCovEffect/"
    UseCovX=c(0,0)
  
  } else if (irun==4){
    ModNm="BetaDevHBMRT_lag1"    
    OutDir="Output/lg1WithCovEffect/"
    UseCovX=c(1,0)
  }
  
  for(itr in 1:2){#1:2
    
    if(RunMainstem==F){
      doTrib=switch(itr,"ubc","lcc")
    } else {
      doTrib=switch(itr,"knights landing","tisdale")
    }
    
    source("GetData.R")

    model_data<-list(Nyrs=Nyrs,Nwks=Nwks,Nestwks=Nestwks,theta=theta,ewk=ewk,Ntot_mu=Ntot_mu,Ntot_sd=Ntot_sd,
                     Nx_mu=Nx_mu,Nx_sd=Nx_sd,CovX=CovX,UseCovX=UseCovX,pr_sd_pro=pr_sd_pro)
                     #Nfor=Nfor,For_ewk=For_ewk,For_Nx_mu=For_Nx_mu,For_Nx_sd=For_Nx_sd,For_CovX=For_CovX,NxError=NxError)
    
    iphi=0.4;ilam=50
    mu_phi=log(iphi/(1-iphi));mu_lambda=log(ilam)
    irt=cbind(rep(mu_phi,Nyrs),rep(mu_lambda,Nyrs))
    
    #inits1=list(muRT=c(mu_phi,mu_lambda),sigmaRT=c(1,1),bCov=c(0,0),sd_pro=rep(0.2,Nyrs),hyp_sd_pro=c(2,5))
    inits1=list(muRT=c(mu_phi,mu_lambda),RTpars=irt,sigmaRT=c(1,1),bCov=c(0,0),sd_pro=1)
    inits2<-inits1;inits3<-inits1;inits<-list(inits1,inits2,inits3)
    
    ParSaveList=c("phi","lambda","cp","sd_pro","RTpars","muRT","sigmaRT","rho","mu_phi","mu_lambda","bCov","For_cp")
    if(irun>=3) ParSaveList=c(ParSaveList,"rho_pro")
    
    Mnm=paste0(ModNm,".stan")
    
    fit=stan(file=Mnm,data=model_data,init=inits,chains=nchains,iter=niter,include=T,pars=ParSaveList)
    
    Snm=paste0(OutDir,doTrib,".Rdata"); save(fit,file=Snm)  
  }#next itr
}#next irun
