
Mnm="JPE"
nchains=3;niter=1000

model_data<-list(Ntribs=Ntribs,Nfor=Nfor, 
                 obs_Nx_mu=obs_Nx_mu, obs_Nx_sd=obs_Nx_sd,
                 srNtot_mu=srNtot_mu, srNtot_sd=srNtot_sd,
                 cp_mu=cp_mu, cp_sd=cp_sd,
                 DS_surv_mu=DS_surv_mu,DS_surv_sd=DS_surv_sd)

#inits1=list(pred_lgNtot=rep(srNtot_mu[itrib],times=Nfor))#,logit_DS_surv=logit(0.1)
inits1=list(pred_lgNtot=matrix(data=rep(srNtot_mu,times=Nfor),nrow=Ntribs,ncol=Nfor))
inits2=list(pred_lgNtot=matrix(data=rep(0.8*srNtot_mu,times=Nfor),nrow=Ntribs,ncol=Nfor))
inits3=list(pred_lgNtot=matrix(data=rep(1.2*srNtot_mu,times=Nfor),nrow=Ntribs,ncol=Nfor))
inits<-list(inits1,inits2,inits3)

#Run stan model and save fit object
Mnm1=paste0(Mnm,".stan")
fit=stan(file=Mnm1,data=model_data,init=inits,chains=nchains,iter=niter,include=T)
#Snm=paste0("Output/Fit_",DoTrib[itrib],".Rdata")
Snm="Output/Fit.Rdata"
save(fit,file=Snm) 

#Save stan model input data for later plotting against poserioes
#Snm=paste0("Output/Data_",DoTrib[itrib],".Rdata")
Snm="Output/Data.Rdata"
save(model_data,file=Snm)

print(summary(fit)$summary)
rhat = summary(fit)$summary[,"Rhat"];ibad=which(rhat[1:length(rhat)-1]>1.05);print("# and vars with rhat>1.05");print(length(ibad));print(rhat[ibad])
print("# of parameters or derived variables saved ");print(length(rhat))

