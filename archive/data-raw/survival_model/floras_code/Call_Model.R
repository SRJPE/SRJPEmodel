library("rstan")
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

#Clear the workspace
rm(list = ls())

source("scripts/GetData_flora.R")

nchains <- 3#3
niter <- 2000

for(irun in 1:9){
  ModNm <- switch(irun,"NoCov","CovWY","CovWY","CovWY","CovWY","CovWY","CovWY","CovWY","CovWY",
                  "CovIndWY","CovIndWY","CovIndWY")
  fnext <- ""
  
  if (ModNm == "NoCov"){
    ParSaveList <- c("P_b","muPb","sdPb","S_bReach","S_bTrib","RE_sd","RE_sdT","S_RE","S_REt","pred_surv",
                     "SurvRelSac","SurvWoodSac","SurvForecast","pred_pcap","pred_survT","TribSurvForecast")
    
    inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
                 S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5,
                 S_bTrib=rep(0,Ntribs),S_REt=rep(-3,NrgT),RE_sdT=0.5)
    
    UseSizeEffect <- 0
    NS_bCov <- 0
    CovX <- rep(0,Nind)
    CovXT <- rep(0,NindT)
    rgwy_ind <- rep(0,Nrg)
    rgwy_indT <- rep(0,NrgT)
    Sz <- rep(0,Nind)
    SzT <- rep(0,NindT)
    Xsz <- rep(0,Nsz)
  
    }else if (ModNm=="CovWY"){ # discrete covariate such as water year type effect + size effect
    ParSaveList <- c("P_b","muPb","sdPb","S_bReach","S_bTrib","S_bCov","S_bCovT","S_bSz","RE_sd","RE_sdT","S_RE","S_REt",
                     "pred_surv","SurvRelSac","SurvWoodSac","SurvForecast","SurvRelSacSz","SurvWoodSacSz","SurvForecastSz",
                     "pred_pcap","pred_survT","pred_survTSz","TribSurvForecast","TribSurvForecastSz")
    
    inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
                   #  S_bCov=rep(0,2), S_bCovT=rep(0,2), #See if removing the initial values altogether works
                   S_bReach=rep(0,3), S_RE=rep(-3,Nrg), RE_sd=0.5)
    
    if(irun==2){ #2 level WY effect and no size effect
    fnext <- "2"
    
    UseSizeEffect <- 0
    CovX <- WY2  
    CovXT <- WY2T
    NS_bCov <- 1
    rgwy_ind <- rgwy2_ind
    rgwy_indT <- rgwy2_indT
    Sz <- rep(0,Nind)
    SzT <- rep(0,NindT)
    Xsz <- rep(0,Nsz)
    
    }else if(irun==3){ #3 level WY effect and no size effect
    fnext <- "3"
    
    UseSizeEffect <- 0
    CovX <- WY3  
    CovXT <- WY3T
    NS_bCov <- 2
    rgwy_ind <- rgwy3_ind
    rgwy_indT <- rgwy3_indT
    Sz <- rep(0,Nind)
    SzT <- rep(0,NindT)
    Xsz <- rep(0,Nsz)
    
    }else if(irun==4){ #2 level WY effect and a fish FL effect
      fnext <- "2_FL"
      
      UseSizeEffect <- 1
      CovX <- WY2  
      CovXT <- WY2T
      NS_bCov <- 1
      rgwy_ind <- rgwy2_ind
      rgwy_indT <- rgwy2_indT
      musz <- mean(c(FL,FL_T),na.rm=TRUE)
      sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
      Sz <- (FL-musz)/sdsz
      SzT <- (FL_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
  
    }else if(irun==5){ #2 level WY effect and a fish weight effect
      fnext <- "2_Wgt"
      
      UseSizeEffect <- 1
      CovX <- WY2  
      CovXT <- WY2T
      NS_bCov <- 1
      rgwy_ind <- rgwy2_ind
      rgwy_indT <- rgwy2_indT
      musz <- mean(c(WGT,WGT_T),na.rm=TRUE)
      sdsz <- sd(c(WGT,WGT_T),na.rm=TRUE)
      Sz <- (WGT-musz)/sdsz
      SzT <- (WGT_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
      
    }else if(irun==6){ #2 level WY effect and a fish condition factor effect
      fnext <- "2_CF"
      
      UseSizeEffect <- 1
      CovX <- WY2  
      CovXT <- WY2T
      NS_bCov <- 1
      rgwy_ind <- rgwy2_ind
      rgwy_indT <- rgwy2_indT
      musz <- mean(c(CF,CF_T),na.rm=TRUE)
      sdsz <- sd(c(CF,CF_T),na.rm=TRUE)
      Sz <- (CF-musz)/sdsz
      SzT <- (CF_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
    
    }else if(irun==7){ #3 level WY effect and a fish FL effect
      fnext <- "3_FL"
      
      UseSizeEffect <- 1
      CovX <- WY3  
      CovXT <- WY3T
      NS_bCov <- 2
      rgwy_ind <- rgwy3_ind
      rgwy_indT <- rgwy3_indT
      musz <- mean(c(FL,FL_T),na.rm=TRUE)
      sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
      Sz <- (FL-musz)/sdsz
      SzT <- (FL_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
      
    }else if(irun==8){ #3 level WY effect and a fish weight effect
      fnext <- "3_Wgt"
      
      UseSizeEffect <- 1
      CovX <- WY3  
      CovXT <- WY3T
      NS_bCov <- 2
      rgwy_ind <- rgwy3_ind
      rgwy_indT <- rgwy3_indT
      musz <- mean(c(WGT,WGT_T),na.rm=TRUE)
      sdsz <- sd(c(WGT,WGT_T),na.rm=TRUE)
      Sz <- (WGT-musz)/sdsz
      SzT <- (WGT_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
      
    }else if(irun==9){ #3 level WY effect and a fish condition factor effect
      fnext <- "3_CF"
      
      UseSizeEffect <- 1
      CovX <- WY3  
      CovXT <- WY3T
      NS_bCov <- 2
      rgwy_ind <- rgwy3_ind
      rgwy_indT <- rgwy3_indT
      musz <- mean(c(CF,CF_T),na.rm=TRUE)
      sdsz <- sd(c(CF,CF_T),na.rm=TRUE)
      Sz <- (CF-musz)/sdsz
      SzT <- (CF_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
      
    }
    
  }else if (ModNm=="CovIndWY"){ # size effect varies with water year type
      ParSaveList <- c("P_b","muPb","sdPb","S_bReach","S_bTrib","S_bCov","S_bCovT","S_bSz","RE_sd","RE_sdT","S_RE","S_REt",
                       "pred_surv","SurvRelSac","SurvWoodSac","SurvForecast","SurvRelSacSz","SurvWoodSacSz","SurvForecastSz",
                       "pred_pcap","pred_survT","pred_survTSz","TribSurvForecast","TribSurvForecastSz")
      
      inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
                     #  S_bCov=rep(0,2), S_bCovT=rep(0,2), #See if removing the initial values altogether works
                     S_bReach=rep(0,3), S_RE=rep(-3,Nrg), RE_sd=0.5)
      
  
      if(irun==10){ #3 level WY effect and a fish length effect that varies with WY
      fnext <- "3_FL"
      
      UseSizeEffect <- 1
      CovX <- WY3  
      CovXT <- WY3T
      NS_bCov <- 2
      rgwy_ind <- rgwy3_ind
      rgwy_indT <- rgwy3_indT
      musz <- mean(c(FL,FL_T),na.rm=TRUE)
      sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
      Sz <- (FL-musz)/sdsz
      SzT <- (FL_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
      
    }else if(irun==11){ #3 level WY effect and a fish weight effect
      fnext <- "3_Wgt"
      
      UseSizeEffect <- 1
      CovX <- WY3  
      CovXT <- WY3T
      NS_bCov <- 2
      rgwy_ind <- rgwy3_ind
      rgwy_indT <- rgwy3_indT
      musz <- mean(c(WGT,WGT_T),na.rm=TRUE)
      sdsz <- sd(c(WGT,WGT_T),na.rm=TRUE)
      Sz <- (WGT-musz)/sdsz
      SzT <- (WGT_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
      
    }else if(irun==12){ #3 level WY effect and a fish condition factor effect
      fnext <- "3_CF"
      
      UseSizeEffect <- 1
      CovX <- WY3  
      CovXT <- WY3T
      NS_bCov <- 2
      rgwy_ind <- rgwy3_ind
      rgwy_indT <- rgwy3_indT
      musz <- mean(c(CF,CF_T),na.rm=TRUE)
      sdsz <- sd(c(CF,CF_T),na.rm=TRUE)
      Sz <- (CF-musz)/sdsz
      SzT <- (CF_T-musz)/sdsz
      Xsz  <- seq(from=min(c(Sz,SzT),na.rm=TRUE),to=max(c(Sz,SzT),na.rm=TRUE),length.out=Nsz)
    }
  }

# else if (ModNm=="CovInd"){#individual covariates such as flow, velocity, temperature + fise size effect
#     ParSaveList=c("P_b","muPb","sdPb","S_bReach","S_bTrib","S_bCov","S_bRT","RE_sd","RE_sdT","S_RE","S_REt","pred_surv",
#                   "SurvWoodSac","SurvWoodSacX","SurvForecast","pred_pcap","pred_survT","TribSurvForecast")
#     
#     inits1<-list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
#                  S_bCov=0, S_bReach=rep(0,Nreaches-1), S_RE=rep(-3,Nrg), RE_sd=0.5)
#     
#     rgwy_ind=rep(0,Nrg);rgwy_indT=rep(0,NrgT) #index for S_bCovT for tributary release groups
#     CovXT=rep(0,NindT)
#     Sz=rep(0,Nind);SzT=rep(0,NindT)
#     
#     if(irun==10){ # flow effect only
#       fnext="_Flow"
#       UseSizeEffect <- 0
#       X <- cbind(d_Sac_sort$Flow_Releasepoint,d_Sac_sort$Flow_Woodson.Br, 
#                  d_Sac_sort$Flow_ButteBridge,d_Sac_sort$Flow_Sacramento)  #Flow covariate at each reach
#      
#       
#     } else if (irun==11){ # velocity effect only
#       fnext="_Vel"
#       UseSizeEffect <- 0
#       
#       X=cbind(d$Vel_Releasepoint,d$Vel_Woodson.Br, d$Vel_ButteBridge,d$Vel_Sacramento) #Velocity
#       
#       
#     } else if (irun==12){ # Temp effect only
#       fnext="_Temp" 
#       UseSizeEffect <- 0
#       
#       X=cbind(d$Temp_Releasepoint,d$Temp_Woodson.Br, d$Temp_ButteBridge,d$Temp_Sacramento) #Temperature
#       
#     } else if(irun==13){ # flow + FL size effects
#       fnext="_Flow_FL"
#       UseSizeEffect <- 1
#       
#       X=cbind(d$Flow_Releasepoint,d$Flow_Woodson.Br, d$Flow_ButteBridge,d$Flow_Sacramento)  #Flow
#       
#       
#     } else if (irun==14){ # temp + FL size effects
#       fnext="_Temp_FL"
#       UseSizeEffect <- 1
#       
#       X=cbind(d$Temp_Releasepoint,d$Temp_Woodson.Br, d$Temp_ButteBridge,d$Temp_Sacramento) #Temperature
#       
#     }
#     
#     Xstd=matrix(nrow=Nind,ncol=Nreaches); muX=mean(X);sdX=sd(X)
#     CovX=matrix(nrow=Nind,ncol=Nreaches)
#     for(j in 1:Nreaches){
#       #muX=mean(X[,j]);sdX=sd(X[,j])#standardize with reach-specific means and standard deviations
#       Xstd[,j]=(X[,j]-muX)/sdX
#       CovX[,j]=Xstd[,j]
#     }
#     Xsz=seq(from=min(Xstd),to=max(Xstd),length.out=Nseq)
# 
# }
 
  model_data<-list(Nind=Nind, Nreaches=Nreaches, Ndetlocs=Ndetlocs,Rmult=Rmult,RmultSac = RmultSac,
                   Nyrs=Nyrs, Nrg=Nrg, CH=CH, yrind=yrind, rgind=rgind,  rgwy_ind=rgwy_ind,rch_covind=rch_covind, 
                   UseSizeEffect = UseSizeEffect,NS_bCov=NS_bCov, CovX=CovX, firstCap=firstCap,lastCap=lastCap,
                   Ntribs=Ntribs,NindT=NindT,NrgT=NrgT,trib_ind=trib_ind,firstCapT=firstCapT,
                   lastCapT=lastCapT,RmultT=RmultT,RmultTrib=RmultTrib,yrindT=yrindT,rgindT=rgindT,CHT=CHT,
                   trib_rg=trib_rg, CovXT=CovXT,rgwy_indT=rgwy_indT,Sz=Sz,SzT=SzT,Nsz=Nsz,Xsz=Xsz)
  
  inits2 <- inits1
  inits3 <- inits1
  inits <- list(inits1,inits2,inits3)
  
  fit <- stan(file=paste0("scripts/",ModNm,".stan"),data=model_data, init=inits, chains=nchains,
              iter=niter,include=T,pars=ParSaveList,seed=1234)#,algorithm = "Fixed_param",sample_file="post.txt",diagnostic_file="diagnostic.txt" #algorithm = "Fixed_param"
  
  fitnm <- paste0("outputs/","fit_",ModNm,fnext,".Rdata")
  
  save(fit,file=fitnm)        

} #next irun


