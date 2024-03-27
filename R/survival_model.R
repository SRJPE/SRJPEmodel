# refactor of Call_Model.R
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

run_survival_model <- function(capture_history_matrix, n_chains = 3, n_iter = 1500, n_years, n_reaches, n_release_groups, n_individuals) {
  # TODO get covariates (in "d" in GetData.R - where is he getting these?)
  # TODO what data are we passing in, and can we get n_years, n_reaches, n_release_groups, n_individuals ?
  # TODO pass in capture_history_matrix, environmental_covariates
  for(irun in 10:12){
    ModNm=switch(irun,"Year_Reach","YR_HBM","NoCov","CovWY2","CovWY2_Reach","CovWY3_Reach","CovInd_Reach","CovInd_Reach","CovInd_Reach","CovInd_ReachSz","CovInd_ReachSz","CovInd_ReachSz")
    fnext=""

    if(ModNm=="Year_Reach"){ #Flora's year*reach model for pCap and survival
      ParSaveList=c("S_b","pred_surv","SurvWoodSac","pred_pcap")
      inits1<-list(S_b=matrix(data=0.85,nrow=Nyrs,ncol=Nreaches),P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches))
      CovX=rep(0,Nind) #no covariate effect so set to 0
      rgwy_ind=rep(0,Nrg)

    } else if (ModNm=="YR_HBM"){
      ParSaveList=c("S_b","P_b","muPb","sdPb","pred_surv","SurvWoodSac","pred_pcap")
      inits1<-list(S_b=matrix(data=0.85,nrow=Nyrs,ncol=Nreaches),P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches))
      CovX=rep(0,Nind) #no covariate effect so set to 0
      rgwy_ind=rep(0,Nrg)

    } else if (ModNm=="NoCov"){
      ParSaveList=c("P_b","muPb","sdPb","S_bReach","RE_sd","S_RE","pred_surv","SurvWoodSac","SurvForecast","pred_pcap")
      inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)
      CovX=rep(0,Nind)
      rgwy_ind=rep(0,Nrg)

    }else if (ModNm=="CovWY2"){ #water year dummy variable (0/1) covariate effect on survival (same for all reaches)
      ParSaveList=c("P_b","muPb","sdPb","S_bReach","S_bCov","RE_sd","S_RE","pred_surv","SurvWoodSac","SurvForecast","pred_pcap")
      inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),S_bCov=1, S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)
      CovX=WY2  #2-level covariate effect (dry or wet)
      rgwy_ind=rep(0,Nrg)

    }else if (ModNm=="CovWY2_Reach"){ #reach-specific water year dummy variable covariate effect on survival (dry vs wet)
      ParSaveList=c("P_b","muPb","sdPb","S_bReach","S_bCov","RE_sd","S_RE","pred_surv","SurvWoodSac","SurvForecast","pred_pcap")
      inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),S_bCov=rep(1,Nreaches-1), S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)
      CovX=WY2 #2 level covariate effect
      rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0) #index for S_bCov for each release group - year (0=dry, 1=wet)

    }else if (ModNm=="CovWY3_Reach"){ #reach-specific 3 water year type covariate effect effect on survival (C, D/BN, W)
      ParSaveList=c("P_b","muPb","sdPb","S_bReach","S_bCov","RE_sd","S_RE","pred_surv","SurvWoodSac","SurvForecast","pred_pcap")
      inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),S_bCov=matrix(data=1,nrow=2,Nreaches-1), S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)
      CovX=WY3  #3 level covariate effect
      rgwy_ind=c(1,1,0,1,2,2,2,1,1,2,1,1,0) #index for S_bCov for each release group - year (0=critical, 1=dry/BN, 2 = wet)

    }else if (ModNm=="CovInd_Reach"){ #Individual covariate effects that vary by reach (flow, temperature, travel time)
      ParSaveList=c("P_b","muPb","sdPb","S_bReach","S_bCov","RE_sd","S_RE","pred_surv","SurvWoodSac","SurvForecast","pred_pcap")
      inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),S_bCov=0, S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)

      if(irun==7){
        X=cbind(d$Flow_Releasepoint,d$Flow_Woodson.Br, d$Flow_ButteBridge,d$Flow_Sacramento)  #Flow
        fnext="_Q"
      } else if (irun==8){
        X=cbind(d$Vel_Releasepoint,d$Vel_Woodson.Br, d$Vel_ButteBridge,d$Vel_Sacramento) #Velocity
        fnext="_Vel"
      } else if (irun==9){
        X=cbind(d$Temp_Releasepoint,d$Temp_Woodson.Br, d$Temp_ButteBridge,d$Temp_Sacramento) #Temperature
        fnext="_Wt"
      }
      Xstd=matrix(nrow=Nind,ncol=Nreaches); muX=mean(X);sdX=sd(X)
      for(j in 1:Nreaches){Xstd[,j]=(X[,j]-muX)/sdX}
      Xseq=seq(from=min(Xstd),to=max(Xstd),length.out=Nseq)

    }else if (ModNm=="CovInd_ReachSz"){ #Individual covariate effects that vary by reach (flow, temperature, travel time)
      ParSaveList=c("P_b","muPb","sdPb","S_bReach","S_bCov","S_bSz","RE_sd","S_RE","pred_surv","SurvWoodSac","SurvWoodSacSz","SurvForecast","SurvForecastSz","pred_pcap")
      inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),S_bCov=0, S_bSz=0, S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)

      if(irun==10){
        X=cbind(d$Vel_Releasepoint,d$Vel_Woodson.Br, d$Vel_ButteBridge,d$Vel_Sacramento) #Velocity
        fnext="_Vel_FL"
        Sz=(FL-mean(FL))/sd(FL)
      } else if (irun==11){
        X=cbind(d$Vel_Releasepoint,d$Vel_Woodson.Br, d$Vel_ButteBridge,d$Vel_Sacramento) #Velocity
        fnext="_Vel_WGT"
        Sz=(WGT-mean(WGT))/sd(WGT)
      } else if (irun==12){
        X=cbind(d$Vel_Releasepoint,d$Vel_Woodson.Br, d$Vel_ButteBridge,d$Vel_Sacramento) #Velocity
        fnext="_Vel_CF"
        Sz=(CF-mean(CF))/sd(CF)
      }

      Xstd=matrix(nrow=Nind,ncol=Nreaches); muX=mean(X);sdX=sd(X)
      for(j in 1:Nreaches){Xstd[,j]=(X[,j]-muX)/sdX}
      Xseq=seq(from=min(Xstd),to=max(Xstd),length.out=Nseq)
      Xsz=seq(from=min(Sz),to=max(Sz),length.out=Nseq)

      CovX=Xstd
      rgwy_ind=rep(0,Nrg)

    }

    if(irun<=6){
      model_data<-list(Nind=Nind, Nreaches=Nreaches, Ndetlocs=Ndetlocs,Rmult=Rmult,Nyrs=Nyrs, Nrg=Nrg, CH=CH, yrind=yrind, rgind=rgind, rch_covind=rch_covind, CovX=CovX, firstCap=firstCap,lastCap=lastCap, rgwy_ind=rgwy_ind)
    } else if (irun<=9) {
      model_data<-list(Nind=Nind, Nreaches=Nreaches, Ndetlocs=Ndetlocs,Rmult=Rmult,Nyrs=Nyrs, Nrg=Nrg, CH=CH, yrind=yrind, rgind=rgind, rch_covind=rch_covind, CovX=CovX, firstCap=firstCap,lastCap=lastCap, rgwy_ind=rgwy_ind,Nseq=Nseq,Xseq=Xseq)
    } else {
      model_data<-list(Nind=Nind, Nreaches=Nreaches, Ndetlocs=Ndetlocs,Rmult=Rmult,Nyrs=Nyrs, Nrg=Nrg, CH=CH, yrind=yrind, rgind=rgind, rch_covind=rch_covind, CovX=CovX, firstCap=firstCap,lastCap=lastCap, rgwy_ind=rgwy_ind,Nseq=Nseq,Xseq=Xseq,Sz=Sz,Xsz=Xsz)
    }


    inits2<-inits1;inits3<-inits1;inits=list(inits1,inits2,inits3)

    fit=stan(file=paste0(ModNm,".stan"),data=model_data, init=inits, chains=nchains,iter=niter,include=T,pars=ParSaveList)#,algorithm = "Fixed_param",sample_file="post.txt",diagnostic_file="diagnostic.txt" #algorithm = "Fixed_param"
    fitnm=paste0(OutDir,"fit_",ModNm,fnext,".Rdata")
    save(fit,file=fitnm)

  } #next irun

  #If fitting model by optimization. Won't work with random effects
  #m <- stan_model(file=paste0(ModNm,".stan"))
  #opt_fit=optimizing(object=m,data=model_data, init=inits1,hessian=T,verbose=T)
  #ipars=which(strtrim(names(opt_fit$par),width=11)=="SurvWoodSac")#example of how to pull out some specific parameter names
  #print(opt_fit$par[ipars])
  #Note the hessian is only available for directly estimated parameters so won't work for SurvWoodSac[] derived variable
  #se=sqrt(diag(solve(opt_fit$hessian/opt_fit$value))); #I think divide by value since optimizer is looking at log probability, not probability. Fuzzy on this bit!
  #print(cbind(opt_fit$par[ipars],se[ipars]))

}


