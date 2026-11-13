#Forecast outmigrant abundance by tributary

Nscale=0.001

print("SR")
Ntribs=length(DoTrib)
predSR_stats=matrix(nrow=Ntribs,ncol=3)
srNtot_mu=vector(length=Ntribs);srNtot_sd=srNtot_mu

for(itrib in 1:Ntribs){
  
  #Get posterior distributions of SR parameters
  Snm=paste0(SRdir,"Fit_",DoTrib[itrib],"_",Am[itrib],"_",SRCovNm[itrib],".Rdata")
  load(file=Snm) 
  dp=as.data.frame(fit, pars = c("alpha","beta","sd_pro"))
  Nsims=dim(dp)[1]
  if(itrib==1)predSR=matrix(data=0,nrow=Nsims,ncol=Ntribs)
  
  icol=which(names(dp)=="alpha");alpha=dp[,icol]; icol=which(names(dp)=="beta");beta=dp[,icol]; icol=which(names(dp)=="sd_pro");sd_pro=dp[,icol]
  
  set.seed(123)#same random # sequence for each tributary
  #srdev=rnorm(n=Nsims,mean=0,sd=mean(sd_pro))#simulate process error for a future single year
  #Best approach which accounts for uncertainty in sd_pro but little difference relative to srdev=rnorm(n=Nsims,mean=0,sd=mean(sd_pro))
  srdev=vector(length=Nsims)
  for(isim in 1:Nsims) srdev[isim]=rnorm(n=1,mean=0,sd=sd_pro[isim])
  
  if(SRCovNm[itrib]=="null"){
    predSR[,itrib]=Sp[itrib]*exp(alpha + beta*Sp[itrib] + srdev)*Nscale
  } else {
    dp=as.data.frame(fit, pars = c("gamma"))
    
    if(SRmodType[itrib]=="Cont"){
      icol=which(names(dp)=="gamma");gamma=dp[,icol]
      predSR[,itrib]=Sp[itrib]*exp(alpha + beta*Sp[itrib] + gamma*Xsr[itrib] + srdev)*Nscale
      
    } else {
      icat=Xsr[itrib]
      if(icat==0){
        gamma=rep(0,times=Nsims)
      } else {
        icol=which(names(dp)==paste0("gamma[",icat,"]"));gamma=dp[,icol]
      }
      predSR[,itrib]=Sp[itrib]*exp(alpha + beta*Sp[itrib] + gamma + srdev)*Nscale
    }
  }
  
  srNtot_mu[itrib]=mean(log(predSR[,itrib]))#mean and sd of of predicted annual abundance for Ntot_Est.stan
  srNtot_sd[itrib]=sd(log(predSR[,itrib]))
  
  predSR_stats[itrib,]=as.double(quantile(predSR[,itrib],probs=CI))#;predSR_stats[itrib,2]=mean(predSR[,itrib])
  print(c(DoTrib[itrib],round(predSR_stats[itrib,],digits=0)))

}#next trib



