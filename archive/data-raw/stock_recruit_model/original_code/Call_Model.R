#MultiRun=F
if(MultiRun==F){
  library("rstan")
  options(mc.cores=parallel::detectCores())
  rstan_options(auto_write=TRUE)
  nchains=3;niter=1000
  
  Mnm="Ricker.stan"#assumess error around lgR/S is normally distributed
  
  OutDir="Output/"
  #DoSite="ubc"
  DoSite="lcc"
  #Am="us"#specify adult data type
  Am="rd"
  #Recommended adult data types by trip from JPE adult report.docx
  #Battle & Clear = upstream passage, redd #Butte = holding, carcass #Deer = holding #Feather = redd or carcass #Mill= redd #Yuba = upstream passage or carcass
  
  #UseCovar=T; UseSQcovar=F #fit covariate model to all years of available covariate data (that overlap with SR data)
  #UseCovar=T; UseSQcovar=T #fit covariate model to years of covariate data only available for all covariates (and that overlap with SR data)
  CoVarNm="rr_max_flow"#appended lifestage and covariate_structure fields. UseCovar=T and UseSQcovar can be T or F
  
  #Fit to years where all covariates are available, but don't use a covariate in fitting (null model for multi-model comparison)
  #CoVarNm=NA; UseCovar=T; UseSQcovar=T
  
  #If wanting to fit model to all brood years and not includding a covariate effect
  #CoVarNm="All";UseCovar=F;UseSQcovar=F; #this will create an SR output file name with All to let you know this is the full broood yr fit
  
}

source("GetData.R")

model_data<-list(Nyrs=Nyrs,SP=SP,mu_obslgRS=mu_obslgRS,sd_obslgRS=sd_obslgRS,X=X)

#Initialization
if(UseCovar==F | (UseCovar==T & is.na(CoVarNm)==T)){
  reg=lm(mu_obslgRS~SP)
  ini_gamma=0
} else {
  reg=lm(mu_obslgRS~SP+X)
  ini_gamma=as.double(reg$coefficients[3])
}
ini_alpha=as.double(reg$coefficients[1]);ini_beta=as.double(reg$coefficients[2])

inits1=list(alpha=ini_alpha,beta=ini_beta,gamma=ini_gamma,sd_pro=1)
inits2=list(alpha=ini_alpha*1.5,beta=ini_beta*0.75,gamma=ini_gamma*0.8,sd_pro=1.25)
inits3=list(alpha=ini_alpha*0.75,beta=ini_beta*1.25,gamma=ini_gamma*1.2,sd_pro=0.75)
inits<-list(inits1,inits2,inits3)

#Call model and save fit object
fit=stan(file=Mnm,data=model_data,init=inits,chains=nchains,iter=niter,include=T)
Snm=paste0(OutDir,"Fit_",DoSite,"_",Am,"_",CoVarNm,".Rdata"); save(fit,file=Snm) 
