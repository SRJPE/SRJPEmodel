library("rstan")


graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

 
RunMainstem=T

#OutDir="Output/NoCovEffect/"
#OutDir="Output/WithCovEffect/"
#parlist=c("muRT[1]","muRT[2]","sigmaRT[1]","sigmaRT[2]","rho","bCov[1]","sd_pro","mu_phi","mu_lambda")
#OutDir="Output/lg1NoCovEffect/"
OutDir="Output/lg1WithCovEffect/"
parlist=c("muRT[1]","muRT[2]","sigmaRT[1]","sigmaRT[2]","rho","bCov[1]","sd_pro","mu_phi","mu_lambda","rho_pro")

print(OutDir)
#parlist=c("mu_phi","mu_lambda","phi","lambda","sd_pro")

npars=length(parlist)
fnstat=paste0(OutDir,"parstats.out")
write(file=fnstat,x=c("site",parlist),ncolumns=npars+1,append=F)
for(itr in 1:2){
  if(RunMainstem==F){
    doTrib=switch(itr,"ubc","lcc")
  } else {
    doTrib=switch(itr,"knights landing","tisdale")
  }
  source("GetData.R")
  fnm=paste0(OutDir,doTrib,".Rdata")
  load(file=fnm)
  print(doTrib)
  print(summary(fit,pars=parlist)$summary)
  xout=itr
  for(j in 1:npars){#median and sd
  	strout=paste0("'",round(summary(fit,pars=parlist)$summary[j,6],digits=2)," (",round(summary(fit,pars=parlist)$summary[j,3],digits=2),")'")
  	xout=c(xout,strout)
  }
  write(file=fnstat,xout,ncolumns=npars+1,append=T)
}

#plot(fit, pars=parlist)
plot(fit, plotfun = "trace", inc_warmup = TRUE,pars=parlist)
#plot(fit, plotfun = "hist

#Check convergence
#niter=fit@stan_args[[1]]$iter;print("# of iterations");print(niter) # total # of iterations
#neff = summary(fit)$summary[,"n_eff"];ibad=which(neff[1:length(neff)-1]<100);print("# and vars with neff<100");print(length(ibad));print(neff[ibad])
#rhat = summary(fit)$summary[,"Rhat"];ibad=which(rhat[1:length(rhat)-1]>1.05);print("# and vars with rhat>1.05");print(length(ibad));print(rhat[ibad])
