library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

OutDir="Output/"
parlist=c("alpha","beta","gamma","sd_pro")

npars=length(parlist)
fnstat=paste0(OutDir,"parstats.out")
write(file=fnstat,x=c("site",parlist),ncolumns=npars+1,append=F)

for(itype in 7:7){
  #DoSite=switch(itype,"ubc","ubc","lcc","lcc")
  #Am=switch(itype,"rd","us","rd","us")
  #CoVarNm=switch(itype,"si_weekly_max_temp_mean","si_max_flow","rr_mean_flow","rr_median_flow")
  #CoVarNm="null"
  
  DoSite=switch(itype,"ubc","ubc","lcc","lcc","mill creek","deer creek","okie dam","knights landing")
  Am=switch(itype,"rd","us","rd","us","rd","hd","hd","us")
  #CoVarNm=switch(itype,"si_max_flow","si_max_flow","rr_mean_flow","si_median_flow","si_min_flow","mi_gdd_sacramento")
  CoVarNm=switch(itype,"si_weekly_max_temp_mean","si_max_flow","rr_mean_flow","rr_median_flow",
                 "si_above_13_temp_week","si_above_13_temp_week","rr_min_flow","rr_max_flow")
  
  Snm=paste0(OutDir,"Fit_",DoSite,"_",Am,"_",CoVarNm,".Rdata")
  load(file=Snm) 
  
  print(c(DoSite,Am,CoVarNm))
  print(summary(fit,pars=parlist)$summary)
  
  xout=Snm
  for(j in 1:npars){#median and sd
    dig=2
    if(j==2) dig=5
    strout=paste0("'",round(summary(fit,pars=parlist)$summary[j,6],digits=dig)," (",round(summary(fit,pars=parlist)$summary[j,3],digits=dig),")'")
    xout=c(xout,strout)
  }
  write(file=fnstat,xout,ncolumns=npars+1,append=T)
}

plot(fit, plotfun = "trace", inc_warmup = TRUE,pars=parlist)
#plot(fit, plotfun = "hist

#Check convergence
#niter=fit@stan_args[[1]]$iter;print("# of iterations");print(niter) # total # of iterations
#neff = summary(fit)$summary[,"n_eff"];ibad=which(neff[1:length(neff)-1]<100);print("# and vars with neff<100");print(length(ibad));print(neff[ibad])
#rhat = summary(fit)$summary[,"Rhat"];ibad=which(rhat[1:length(rhat)-1]>1.05);print("# and vars with rhat>1.05");print(length(ibad));print(rhat[ibad])
