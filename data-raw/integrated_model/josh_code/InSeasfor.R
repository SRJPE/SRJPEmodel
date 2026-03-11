predIn_stats=array(dim=c(Ntribs,Nfor,3))
cp_mu=matrix(nrow=Ntribs,ncol=Nfor);cp_sd=cp_mu
obs_Nx_mu=cp_mu;obs_Nx_cv=cp_mu;obs_Nx_sd=cp_mu

print("Inseason")
for(itrib in 1:Ntribs){
  for(ifor in 1:Nfor){
    
    For_ewk=which(calwk==For_CalWk[ifor])
    
    #Currently setup to use null model based on proportion of run forecasts caclculated in stan model.
    #Will eventually replace this with code to use covariate model specified in For_Input.csv
    #This will require reading in parameters and calculating For_cp within this code
    fnm=paste0(InseasDir,DoTrib[itrib],"_null.Rdata")
    load(file=fnm)
    
    dis=as.data.frame(fit, pars = c("For_cp"))
    Nsims=dim(dis)[1]
    
    cp_all=dis[,For_ewk];grecs=which(is.finite(cp_all)==T & cp_all>0 & cp_all<1)
    Ngsims=length(grecs)
    cp=cp_all[grecs]
    
    cp_mu[itrib,ifor]=mean(logit(cp));cp_sd[itrib,ifor]=sd(logit(cp))

    obs_Nx_mu[itrib,ifor]=log(For_Nx_mu[itrib,ifor])#total abundance through this week in log space
    obs_Nx_cv[itrib,ifor]=For_Nx_sd[itrib,ifor]/For_Nx_mu[itrib,ifor]#the cv in untranformed space
    obs_Nx_sd[itrib,ifor]=sqrt(log(obs_Nx_cv[itrib,ifor]^2+1))#convert from cv in untransformed space to sd in log space
    
    #Uncertainty in cummulative abundance on forecast week in log space
    pNx=rnorm(n=Ngsims,mean=obs_Nx_mu[itrib,ifor],sd=obs_Nx_sd[itrib,ifor])
    
    #estimate of annual abundance (cum abundance on forecast wk expaned by proportion of run by that wk)
    predIn=exp(pNx)/cp 
    #predIn=mean(exp(pNx))/cp #no uncertainty
    
    predIn_stats[itrib,ifor,]=as.double(quantile(predIn,probs=CI))#;predIn_stats[itrib,ifor,2]=mean(predIn)
    print(c(DoTrib[itrib],For_ewk,round(predIn_stats[itrib,ifor,],digits=0)))
  }#next ifor
}

