library("rstan")
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

#Clear the workspace
rm(list = ls())

source("scripts/GetData.R")

nchains <- 3 #3
niter <- 1000
for(irun in 26:26){
  ModNm <- switch(irun,"NoCov","CovIndWY","CovIndWY","CovIndWY","CovIndWY",
                  "CovIndWY","CovIndWY","CovIndWY","CovIndWY","CovIndWY",
                  "CovIndWY","CovIndWY","CovIndWY","CovIndWY","CovIndCont",
                  "CovIndCont","CovIndCont","CovIndCont","CovIndCont","CovIndCont")
  fnext <- ""

  if (ModNm == "NoCov"){
    ParSaveList <- c("P_b","muPb","sdPb","S_bReach","S_bTrib","RE_sd","RE_sdT","S_RE","S_REt","log_lik",
                     "pred_surv","SurvRelSac","SurvForecast","pred_pcap","pred_survT","TribSurvForecast",
                     "SurvForecast_nore","TribSurvForecast_nore")

    inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),
                   muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches))
                 #S_bReach=rep(0,Nreaches), S_RE=rep(-3,Nrg), RE_sd=0.5,
                 #S_bTrib=rep(0,Ntribs),S_REt=rep(-3,NrgT),RE_sdT=0.5)

    UseSizeEffect <- 0
    NS_bCov <- 0
    CovX <- rep(0,Nind)
    CovXT <- rep(0,NindT)
    rgwy_ind <- rep(0,Nrg)
    rgwy_indT <- rep(0,NrgT)
    Sz <- rep(0,Nind)
    SzT <- rep(0,NindT)
    Xsz <- rep(0,Nsz)
    Xvec <- rep(0,NsX)
    XvecT <- rep(0,NsX)

    }else if (ModNm=="CovIndWY_0820"){ #discrete variable such as water year type effect + a fish condition variable
      ParSaveList <- c("TT_bCov","TT_bCovT","TT_RE","TT_RET","TTRE_sd","TTRE_sdT",
                       "Pro_sd","Pro_sdT","T_bSz","TT_reach","TT_RelSac","TT_reachT",
                       "P_b","muPb","sdPb","S_bReach","S_bTrib","T_bReach","T_bTrib","S_bCov","S_bCovT","S_bSz",
                       "RE_sd","RE_sdT","S_RE","S_REt",'log_lik',
                       "pred_surv","SurvRelSac","SurvRelSacSz","pred_survT","pred_survTSz",
                       "pred_surv_per100","pred_survT_per100","pred_pcap",
                       "SurvForecast","SurvForecastSz","SurvForecast_nore","SurvForecastSz_nore",
                       "TribSurvForecast","TribSurvForecastSz","TribSurvForecast_nore","TribSurvForecastSz_nore",
                       "SurvForecastSz_rst","TribSurvForecastSz_rst","TTForecastSz_rst","TTForecastTSz_rst",
                       "TTForecast","TTForecastT","TTForecastSz","TTForecastTSz")

    inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches))
                   #  S_bCov=rep(0,2), S_bCovT=rep(0,2), #See if removing the initial values altogether works
                   #S_bReach=rep(0,Nreaches), S_RE=rep(-3,Nrg), RE_sd=0.5)

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
    Xvec <- rep(0,NsX)
    XvecT <- rep(0,NsX)

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
    Xvec <- rep(0,NsX)
    XvecT <- rep(0,NsX)

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
      Xsz  <-  (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

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
      Xsz  <- (seq(from=4,to=35,length.out=Nsz)-musz)/sdsz
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

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
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

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
      Xsz  <-  (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

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
      Xsz  <- (seq(from=4,to=35,length.out=Nsz)-musz)/sdsz
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

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
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

    }else if(irun==10){ #No WY effect and a fish length effect
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
      Xsz  <- (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

    }else if(irun==11){ #No WY effect and a fish weight effect
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
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

    }else if(irun==12){ #No WY effect and a fish condition factor effect
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
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

    }else if(irun==13){ #
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
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)

    }else if(irun==14){ #
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
      Xsz  <-  (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
      Xvec <- rep(0,NsX)
      XvecT <- rep(0,NsX)
    }

    }else if (ModNm=="CovIndCont"){ # continuous variable at the individual fish level such as flow at release + a fish condition variable
      ParSaveList <- c("TT_bCov","TT_bCovT","TT_RE","TT_RET","TTRE_sd","TTRE_sdT",
                       "Pro_sd","Pro_sdT","T_bSz","TT_reach","TT_RelSac","TT_reachT",
                      "P_b","muPb","sdPb","S_bReach","S_bTrib","T_bReach","T_bTrib","S_bCov","S_bCovT","S_bSz",
                      "RE_sd","RE_sdT","S_RE","S_REt",'log_lik',
                      "pred_surv","SurvRelSac","SurvRelSacSz","pred_survT","pred_survTSz",
                      "pred_surv_per100","pred_survT_per100","pred_pcap",
                      "SurvForecast","SurvForecastSz","SurvForecast_nore","SurvForecastSz_nore",
                      "TribSurvForecast","TribSurvForecastSz","TribSurvForecast_nore","TribSurvForecastSz_nore",
                      "SurvForecastSz_rst","TribSurvForecastSz_rst","TTForecastSz_rst","TTForecastTSz_rst",
                      "TTForecast","TTForecastT","TTForecastSz","TTForecastTSz")

      inits1 <- list(P_b=matrix(data=2.2,nrow=Nyrs,ncol=Nreaches),muPb=rep(0,Nreaches),sdPb=rep(1.0,Nreaches),
                     #  S_bCov=rep(0,2), S_bCovT=rep(0,2), S_RE=rep(-3,Nrg), RE_sd=0.5 #See if removing the initial values altogether works
                     S_bReach=c(1.5,0.5,0.5,0.5) )

      if (irun==15){ #Flow at RBDD for each sac fish and at release for Butte and Feather fish
        fnext <- "_MaxFlow_FL"

        UseSizeEffect <- 1
        NS_bCov <- 0
        CovX <- data.frame(cbind(MaxflowSac.z,MaxflowSac.z,MaxflowSac.z,MaxflowDelta.z))
        CovXT <-data.frame(cbind(MaxflowT.z,MaxflowDeltaT.z))
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
        Xsz  <- (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
        Xvec <- (seq(from =4000,to=40000,length.out=NsX)-mux)/sdx
        XvecT <- (seq(from =150,to=30500,length.out=NsX)-muxT)/sdxT

      } else if (irun==16){ #Flow at RBDD for each sac fish and at release for Butte and Feather fish
        fnext <- "_MaxFlow"

        UseSizeEffect <- 0
        NS_bCov <- 0
        CovX <- data.frame(cbind(MaxflowSac.z,MaxflowSac.z,MaxflowSac.z,MaxflowDelta.z))
        CovXT <-data.frame(cbind(MaxflowT.z,MaxflowDeltaT.z))
        mux <- mean(MaxflowSac,na.rm=TRUE)
        sdx <- sd(MaxflowSac,na.rm=TRUE)
        muxT <- mean(MaxflowT,na.rm=TRUE)
        sdxT <- sd(MaxflowT,na.rm=TRUE)
        rgwy_ind <- rep(0,Nrg)
        rgwy_indT <- rep(0,NrgT)
        Sz <- rep(0,Nind)
        SzT <- rep(0,NindT)
        Xsz <- rep(0,Nsz)
        Xvec <- (seq(from =4000,to=40000,length.out=NsX)-mux)/sdx
        XvecT <- (seq(from =150,to=30500,length.out=NsX)-muxT)/sdxT

      } else if (irun==17) { #Temp at RBDD for each sac fish and at release for Butte and Feather fish
        fnext <- "_TempRel"

        UseSizeEffect <- 0
        NS_bCov <- 0
        CovX <- TempRel
        CovXT <-  TempRelT
        rgwy_ind <- rep(0,Nrg)
        rgwy_indT <- rep(0,NrgT)
        Sz <- rep(0,Nind)
        SzT <- rep(0,NindT)
        Xsz <- rep(0,Nsz)
        Xvec <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)
        XvecT <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)

      } else if (irun==18){#Temp at RBDD for each sac fish and at release for Butte and Feather fish + FL effect
        fnext <- "_TempRel_FL"

        UseSizeEffect <- 1
        NS_bCov <- 0
        CovX <- TempRel
        CovXT <-  TempRelT
        rgwy_ind <- rep(0,Nrg)
        rgwy_indT <- rep(0,NrgT)
        musz <- mean(c(FL,FL_T),na.rm=TRUE)
        sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
        Sz <- (FL-musz)/sdsz
        SzT <- (FL_T-musz)/sdsz
        Xsz  <- seq(from=30,to=150,length.out=Nsz)
        Xvec <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)
        XvecT <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)

      } else if (irun==19){ # Add date of release fixed effect term
        fnext <- "_DoR_FL"

        UseSizeEffect <- 1
        NS_bCov <- 0
        CovX <- DoR
        CovXT <-  DoRT
        mux <- mean(c(DoR,DoRT),na.rm=TRUE)
        sdx <- sd(c(DoR,DoRT),na.rm=TRUE)
        rgwy_ind <- rep(0,Nrg)
        rgwy_indT <- rep(0,NrgT)
        musz <- mean(c(FL,FL_T),na.rm=TRUE)
        sdsz <- sd(c(FL,FL_T),na.rm=TRUE)
        Sz <- (FL-musz)/sdsz
        SzT <- (FL_T-musz)/sdsz
        Xsz  <- (seq(from=10,to=150,length.out=Nsz)-musz)/sdsz
        Xvec <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)
        XvecT <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)

    } else if (irun==20){ # Add date of release fixed effect term
      fnext <- "_DoR"

      UseSizeEffect <- 0
      NS_bCov <- 0
      CovX <- DoR
      CovXT <-  DoRT
      mux <- mean(c(DoR,DoRT),na.rm=TRUE)
      sdx <- sd(c(DoR,DoRT),na.rm=TRUE)
      rgwy_ind <- rep(0,Nrg)
      rgwy_indT <- rep(0,NrgT)
      Sz <- rep(0,Nind)
      SzT <- rep(0,NindT)
      Xsz <- rep(0,Nsz)
      Xvec <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)
      XvecT <- seq(from=min(c(CovX,CovXT),na.rm=TRUE),to=max(c(CovX,CovXT),na.rm=TRUE),length.out=NsX)
    }
    }

  model_data<-list(Nind=Nind,Nreaches=Nreaches,Ndetlocs=Ndetlocs,Rmult=Rmult,RmultSac=RmultSac,RmultTis=RmultTis,
                   Rmultrst=Rmultrst,Nrst=Nrst,RmultTrst=RmultTrst,NrstT=NrstT,TribForRST=TribForRST,
                   ReachKMrst=ReachKMrst,ReachKMTrst=ReachKMTrst,
                   Nyrs=Nyrs,Nrg=Nrg,CH=CH,yrind=yrind, rgind=rgind,rgwy_ind=rgwy_ind,rch_covind=rch_covind,
                   UseSizeEffect=UseSizeEffect,NS_bCov=NS_bCov,CovX=CovX,firstCap=firstCap,lastCap=lastCap,
                   Ntribs=Ntribs,NindT=NindT,NreachesT=NreachesT,NrgT=NrgT,trib_ind=trib_ind,firstCapT=firstCapT,
                   lastCapT=lastCapT,RmultT=RmultT,RmultTrib=RmultTrib,yrindT=yrindT,rgindT=rgindT,CHT=CHT,
                   trib_rg=trib_rg,CovXT=CovXT,rgwy_indT=rgwy_indT,Sz=Sz,SzT=SzT,Nsz=Nsz,
                   Xsz=Xsz,NsX=NsX,Xvec=Xvec,XvecT=XvecT,
                   Nobs=Nobs,NobsT=NobsT, ObsTT=ObsTT,TTind=TTind,ObsTTT=ObsTTT, TTindT=TTindT,
                  ReachKM=ReachKM,ReachKMT=ReachKMT,ReachKM_ind=ReachKM_ind,ReachKMT_ind=ReachKMT_ind)

  inits2 <- inits1
  inits3 <- inits1
  inits <- list(inits1,inits2,inits3)

  fit <- stan(file=paste0("scripts/",ModNm,".stan"),data=model_data, init=inits, chains=nchains,
              iter=niter,include=T,pars=ParSaveList,seed=1234)#,algorithm = "Fixed_param",sample_file="post.txt",diagnostic_file="diagnostic.txt" #algorithm = "Fixed_param"

  fitnm <- paste0("outputs/","fit_",ModNm,fnext,".Rdata")

  save(fit,file=fitnm)

} #next irun


