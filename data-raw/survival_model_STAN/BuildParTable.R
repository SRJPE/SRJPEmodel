#The program dumps some critical parameters values to console to be copied into CovSum sheet in Summary.xlsx

library("rstan")

varnm1="S_bCov";varnm2="RE_sd"

for(irun in 1:6){
  ModNm=switch(irun,"Year_Reach","CovWY2","CovWY2_Reach","CovWY3_Reach","CovInd_Reach","CovInd_Reach","CovInd_Reach")

  fnext=""
  if(irun==5){
    fnext="_Q"
  } else if (irun==6){
    fnext="_Vel"
  } else if (irun==7){
    fnext="_Wt"
  }
  
  load(file=paste0("Results/fit_",ModNm,fnext,".Rdata"))
  
  mu_S_bReach=vector(length=3);sd_S_bReach=mu_S_bReach
  S_bReach=as.data.frame(fit,pars=c("S_bReach"))
  for(j in 1:4){
    vname=paste0("S_bReach[",j,"]")
    icol=which(names(S_bReach)==vname)
    mu_S_bReach[j]=round(mean(S_bReach[,icol]),digits=2)
    sd_S_bReach[j]=round(sd(S_bReach[,icol]),digits=2)
    
  }

  S_bCov=as.data.frame(fit, pars = c(varnm1))  
  if(irun==2 | irun==5 | irun==6 | irun==7){
    mu_S_bCov=round(mean(unlist(S_bCov)),digits=2);sd_S_bCov=round(sd(unlist(S_bCov)),digits=2)
 
  } else if (irun==3){
    mu_S_bCov=vector(length=3);sd_S_bCov=mu_S_bCov
    for(j in 1:3){
      vname=paste0(varnm1,"[",j,"]")
      icol=which(names(S_bCov)==vname)
      mu_S_bCov[j]=round(mean(S_bCov[,icol]),digits=2)
      sd_S_bCov[j]=round(sd(S_bCov[,icol]),digits=2)
    }

  } else if (irun==4){
    mu_S_bCov=matrix(nrow=2,ncol=3);sd_S_bCov=mu_S_bCov
    for(iwy in 1:2){
      for(j in 1:3){
        vname=paste0(varnm1,"[",iwy,",",j,"]")
        icol=which(names(S_bCov)==vname)
        mu_S_bCov[iwy,j]=round(mean(S_bCov[,icol]),digits=2)
        sd_S_bCov[iwy,j]=round(sd(S_bCov[,icol]),digits=2)
      }
    }
  
  }
  
  RE_sd=unlist(as.data.frame(fit, pars = c(varnm2)))
  mu_RE_sd=round(mean(RE_sd),digits=2); sd_RE_sd=round(sd(RE_sd),digits=2)
  
  print("")
  print(ModNm)
  print(mu_S_bReach);print(sd_S_bReach)
  print(mu_S_bCov);print(sd_S_bCov)
  print(mu_RE_sd);print(sd_RE_sd)
} #next irun
