library("rstan")

#load(file="Results/fit_Year_Reach.Rdata")
#load(file="Results/fit_NoCov.Rdata")
#load(file="Results/fit_CovInd_Reach_Q.Rdata")
#load(file="Results/fit_CovWY2.Rdata")
#load(file="Results/fit_CovWY2_Reach.Rdata")
#load(file="Results/fit_CovWY3_Reach.Rdata")
#parlist=c("muPb","sdPb","S_bReach","RE_sd","S_bCov")#,c("S_bReach","S_bCov")#pred_surv","SurvWoodSac")#,"SurvForecast","m_RE","RE_sd""pred_surv","M_bCov",

#load(file="Results/fit_CovInd_Reach_Vel.Rdata")
load(file="Results/fit_CovInd_ReachSz_Vel_Wgt.Rdata")
parlist=c("S_bCov","RE_sd") #,"S_bSz",
print(summary(fit,pars=parlist)$summary)
plot(fit,pars=parlist)
plot(fit,plotfun="hist",pars=parlist)
plot(fit, plotfun = "trace", inc_warmup = TRUE,pars=parlist)
#plot(fit_chlm,plotfun="rhat")

#write.table(file="temp.out",summary(fit,pars=parlist)$summary,row.names=T,col.names=T)

