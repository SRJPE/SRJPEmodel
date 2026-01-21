library("rstan")
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

#Clear the workspace
rm(list = ls())

source("scripts/GetData.R")

nchains <- 3
niter <- 1000
for(irun in 15:15){
  ModNm <- switch(irun,"NoCov","CovIndWY","CovIndWY","CovIndWY","CovIndWY",
                  "CovIndWY","CovIndWY","CovIndWY","CovIndWY","CovIndWY",
                  "CovIndWY","CovIndWY","CovIndWY","CovIndWY",
                  "CovIndCont", ,"CovIndCont","CovIndCont","CovIndCont",
                  "CovIndCont","CovIndCont")
  fnext <- ""
  
  if (ModNm == "NoCov"){
    ParSaveList <- c("P_b","muPb","sdPb","S_bReach","S_bTrib","RE_sd","RE_sdT","S_RE","S_REt","log_lik", 
                     "pred_surv","SurvRelSac","SurvForecast","pred_pcap","pred_survT","TribSurvForecast",
                     "SurvForecast_nore","TribSurvForecast_nore")
    
    inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches))

    UseSizeEffect <- 0
    NS_bCov <- 0
    CovX <- rep(0,Nind)
    CovXT <- rep(0,NindT)
    rgwy_ind <- rep(0,Nrg)
    rgwy_indT <- rep(0,NrgT)
    Sz <- rep(0,Nind)
    SzT <- rep(0,NindT)
    Xsz <- rep(0,Nsz)
    X <- rep(0,Nind)
    XT <- rep(0,NindT)
    Xvec <- rep(0,NsX)
    XvecT <- rep(0,NsX)
  
    }else if (ModNm=="CovIndWY"){ #Individual fish condition variable + discrete variable such as water year (WY) type effect 
    ParSaveList <- c("P_b","muPb","sdPb","S_bReach","S_bTrib","S_bCov","S_bCovT","S_bSz",
                     "RE_sd","RE_sdT","S_RE","S_REt",'log_lik',
                     "pred_surv","SurvRelSac","pred_survT","pred_survTSz","pred_pcap",
                     "SurvForecast","SurvForecastSz","SurvForecast_nore","SurvForecastSz_nore",
                     "TribSurvForecast","TribSurvForecastSz","TribSurvForecast_nore","TribSurvForecastSz_nore")
    
    inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches))
    
    if(irun==2){ #2 level WY and no size factor
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
    X <- rep(0,Nrg)
    XT <- rep(0,NrgT)
    Xvec <- rep(0,NsX)
    XvecT <- rep(0,NsX)
    
    }else if(irun==3){ #3 level WY and no size factor
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
    X <- rep(0,Nrg)
    XT <- rep(0,NrgT)
    Xvec <- rep(0,NsX)
    XvecT <- rep(0,NsX)
    
    }else if(irun==4){ #2 level WY and a fish FL factor
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
      Xsz  <-  (seq(from=10,to=130,length.out=Nsz)-musz)/sdsz
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
  
    }else if(irun==5){ #2 level WY and a fish weight factor
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
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      
    }else if(irun==6){ #2 level WY and a fish condition factor
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
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
    
    }else if(irun==7){ #3 level WY effect and a fish FL factor
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
      Xsz  <-  (seq(from=10,to=130,length.out=Nsz)-musz)/sdsz
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      
    }else if(irun==8){ #3 level WY and a fish weight factor
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
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      
    }else if(irun==9){ #3 level WY and a fish condition factor
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
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      
    }else if(irun==10){ #No WY and a fish length factor
      fnext <- "_FL"
      
      UseSizeEffect <- 1
      NS_bCov <- 0
      CovX <- rep(0,Nind)
      CovXT <- rep(0,NindT)
      rgwy_ind <- rep(0,Nrg)
      rgwy_indT <- rep(0,NrgT)
      musz <- mean(c(FL,FL_T),na.rm=TRUE)
      sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
      Sz <- (FL-musz)/sdsz
      SzT <- (FL_T-musz)/sdsz
      Xsz  <- (seq(from=10,to=130,length.out=Nsz)-musz)/sdsz
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)
      
    }else if(irun==11){ #No WY and a fish weight factor
      fnext <- "_Wgt"
      
      UseSizeEffect <- 1
      NS_bCov <- 0
      CovX <- rep(0,Nind)
      CovXT <- rep(0,NindT)
      rgwy_ind <- rep(0,Nrg)
      rgwy_indT <- rep(0,NrgT)
      musz <- mean(c(WGT,WGT_T),na.rm=TRUE)
      sdsz <- sd(c(WGT,WGT_T),na.rm=TRUE)
      Sz <- (WGT-musz)/sdsz
      SzT <- (WGT_T-musz)/sdsz
      Xsz  <- (seq(from=4,to=35,length.out=Nsz)-musz)/sdsz
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)
      
    }else if(irun==12){ #No WY and a fish condition factor 
      fnext <- "_CF"
      
      UseSizeEffect <- 1
      NS_bCov <- 0
      CovX <- rep(0,Nind)
      CovXT <- rep(0,NindT)
      rgwy_ind <- rep(0,Nrg)
      rgwy_indT <- rep(0,NrgT)
      musz <- mean(c(CF,CF_T),na.rm=TRUE)
      sdsz <- sd(c(CF,CF_T),na.rm=TRUE)
      Sz <- (CF-musz)/sdsz
      SzT <- (CF_T-musz)/sdsz
      Xsz  <- (seq(from= 0.0006,to=0.002,length.out=Nsz)-musz)/sdsz
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)
      
    }else if(irun==13){ #Flow exceedance index from CDFW and no fish factor
      fnext <- "_Fexceed"
      
      UseSizeEffect <- 0
      CovX <- FlowexceedSac 
      CovXT <- FlowexceedT
      NS_bCov <- 2
      rgwy_ind <- flowexceedSac_ind 
      rgwy_indT <- flowexceedT_ind 
      Sz <- rep(0,Nind)
      SzT <- rep(0,NindT)
      Xsz <- rep(0,Nsz)
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)
      
    }else if(irun==14){ #Flow exceedance index and a fish FL factor
      fnext <- "_Fexceed_FL"
      
      UseSizeEffect <- 1
      CovX <- FlowexceedSac 
      CovXT <- FlowexceedT
      NS_bCov <- 2
      rgwy_ind <- flowexceedSac_ind 
      rgwy_indT <- flowexceedT_ind 
      musz <- mean(c(FL,FL_T),na.rm=TRUE)
      sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
      Sz <- (FL-musz)/sdsz
      SzT <- (FL_T-musz)/sdsz
      Xsz  <-  (seq(from=10,to=130,length.out=Nsz)-musz)/sdsz
      X <- rep(0,Nrg)
      XT <- rep(0,NrgT)
      Xvec <- rep(0,NsX)
    
    }else if (ModNm=="CovIndContTT"){ #Individual fish condition variable + continuous variable 
      ParSaveList <- c("TT_bCov","TT_bCovT","lgB0","lgB0T","TT_RE","TT_RET","TTRE_sd","TTRE_sdT",
                       "Pro_sd","Pro_sdT","T_bSz","TT_reach","TT_RelSac","TT_reachT",
                      "P_b","muPb","sdPb","S_bReach","S_bTrib","S_bCov","S_bCovT","S_bSz",
                      "RE_sd","RE_sdT","S_RE","S_REt",'log_lik',
                      "pred_surv","SurvRelSac","SurvRelSacSz","pred_survT","pred_survTSz",
                      "pred_surv_per100","pred_survT_per100","pred_pcap",
                      "SurvForecast","SurvForecastSz","SurvForecast_nore","SurvForecastSz_nore",
                      "TribSurvForecast","TribSurvForecastSz","TribSurvForecast_nore","TribSurvForecastSz_nore",
                      "SurvForecastSz_rst","TribSurvForecastSz_rst","TTForecastSz_rst","TTForecastTSz_rst",
                      "TTForecast","TTForecastT","TTForecastSz","TTForecastTSz")
      
      inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
                    S_bReach=c(1.5,0.5,0.5,0.5) )
        
      if (irun==15){ #Peak flow during the month of release and no fish factor
      fnext <- "_MaxFlow"
      
      UseSizeEffect <- 0
      NS_bCov <- 0
      CovX <- MaxflowSac.z
      CovXT <-  MaxflowT.z
      mux <- mean(MaxflowSac,na.rm=TRUE)
      sdx <- sd(MaxflowSac,na.rm=TRUE)
      muxT <- mean(MaxflowT,na.rm=TRUE)
      sdxT <- sd(MaxflowT,na.rm=TRUE)
      rgwy_ind <- rep(0,Nrg)
      rgwy_indT <- rep(0,NrgT)
      Sz <- rep(0,Nind)
      SzT <- rep(0,NindT)
      Xsz <- rep(0,Nsz)
      X <- maxflowSac.df$ind
      XT <- maxflowT.df$ind
      Xvec <- (seq(from =4000,to=40000,length.out=NsX)-mux)/sdx
      XvecT <- (seq(from =150,to=30500,length.out=NsX)-muxT)/sdxT
        
      } else if (irun==17){ #Peak flow during the month of release and fish FL factor
          fnext <- "_MaxFlow_FL"
          
          UseSizeEffect <- 1
          NS_bCov <- 0
          CovX <- MaxflowSac.z
          CovXT <-  MaxflowT.z
          mux <- mean(MaxflowSac,na.rm=TRUE)
          sdx <- sd(MaxflowSac,na.rm=TRUE)
          muxT <- mean(MaxflowT,na.rm=TRUE)
          sdxT <- sd(MaxflowT,na.rm=TRUE)
          rgwy_ind <- rep(0,Nrg)
          rgwy_indT <- rep(0,NrgT)
          musz <- mean(c(FL,FL_T),na.rm=TRUE)
          sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
          Sz <- (FL-musz)/sdsz
          SzT <- (FL_T-musz)/sdsz
          Xsz  <- (seq(from=10,to=130,length.out=Nsz)-musz)/sdsz
          X <- maxflowSac.df$ind
          XT <- maxflowT.df$ind
          Xvec <- (seq(from =4000,to=40000,length.out=NsX)-mux)/sdx
          XvecT <- (seq(from =150,to=30500,length.out=NsX)-muxT)/sdxT
       }
    }
    
   model_data<-list(Nind=Nind,Nreaches=Nreaches,Ndetlocs=Ndetlocs,Rmult=Rmult,RmultSac=RmultSac,RmultTis=RmultTis,
                   Rmultrst=Rmultrst,Nrst=Nrst,RmultTrst=RmultTrst,NrstT=NrstT,ReachKMrst=ReachKMrst,ReachKMTrst=ReachKMTrst,
                   Nyrs=Nyrs,Nrg=Nrg,CH=CH,yrind=yrind, rgind=rgind,rgwy_ind=rgwy_ind,rch_covind=rch_covind, 
                   UseSizeEffect=UseSizeEffect,NS_bCov=NS_bCov,CovX=CovX,firstCap=firstCap,lastCap=lastCap,
                   Ntribs=Ntribs,NindT=NindT,NreachesT=NreachesT,NrgT=NrgT,trib_ind=trib_ind,firstCapT=firstCapT,
                   lastCapT=lastCapT,RmultT=RmultT,RmultTrib=RmultTrib,yrindT=yrindT,rgindT=rgindT,CHT=CHT,
                   trib_rg=trib_rg,CovXT=CovXT,rgwy_indT=rgwy_indT,Sz=Sz,SzT=SzT,Nsz=Nsz,
                   Xsz=Xsz,NsX=NsX,X=X,XT=XT,Xvec=Xvec,XvecT=XvecT,
                   Nobs=Nobs, NobsT=NobsT, ObsTT=ObsTT,TTind=TTind,ObsTTT=ObsTTT, TTindT=TTindT,
                  ReachKM=ReachKM,ReachKMT=ReachKMT,ReachKM_ind=ReachKM_ind,ReachKMT_ind=ReachKMT_ind)
  
  inits2 <- inits1
  inits3 <- inits1
  inits <- list(inits1,inits2,inits3)
  
  fit <- stan(file=paste0("scripts/",ModNm,".stan"),data=model_data, init=inits, chains=nchains,
              iter=niter,include=T,pars=ParSaveList,seed=1234)#,algorithm = "Fixed_param",sample_file="post.txt",diagnostic_file="diagnostic.txt" #algorithm = "Fixed_param"
  
  fitnm <- paste0("outputs/","fit_",ModNm,fnext,".Rdata")
  
  save(fit,file=fitnm)        

} #next irun


