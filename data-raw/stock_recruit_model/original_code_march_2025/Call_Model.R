
source("GetData.R")

model_data<-list(Nyrs=Nyrs,SP=SP,mu_obslgRS=mu_obslgRS,sd_obslgRS=sd_obslgRS,X=X)

#Initialization
if(CoVarNm=="null" | CoVarNm=="WY2" | CoVarNm=="WY3"){
  reg=lm(mu_obslgRS~SP)
  ini_gamma=0
} else if (CoVarNm!="null" & CoVarNm!="mu" & CoVarNm!="WY2" & CoVarNm!="WY3") {
  reg=lm(mu_obslgRS~SP+X)
  ini_gamma=as.double(reg$coefficients[3])
}

if(CoVarNm=="mu"){
  Mnm="Ricker_mu.stan"
  ini_alpha=mean(mu_obslgRS)
  inits1=list(alpha=ini_alpha,sd_pro=1)
  inits2=list(alpha=ini_alpha*1.5,sd_pro=1)
  inits3=list(alpha=ini_alpha*0.75,sd_pro=1)
  model_data<-list(Nyrs=Nyrs,SP=SP,mu_obslgRS=mu_obslgRS,sd_obslgRS=sd_obslgRS)

} else if(CoVarNm=="null") {
  Mnm="Ricker_null.stan"#assumess error around lgR/S is normally distributed
  ini_alpha=as.double(reg$coefficients[1]);ini_beta=as.double(reg$coefficients[2])
  if(ini_beta>0) ini_beta=-0.0001
  
  inits1=list(alpha=ini_alpha,beta=ini_beta,sd_pro=1)
  inits2=list(alpha=ini_alpha*1.5,beta=ini_beta*0.75,sd_pro=1.25)
  inits3=list(alpha=ini_alpha*0.75,beta=ini_beta*1.25,sd_pro=0.75)
  model_data<-list(Nyrs=Nyrs,SP=SP,mu_obslgRS=mu_obslgRS,sd_obslgRS=sd_obslgRS)

} else if (CoVarNm=="WY2" | CoVarNm=="WY3"){
  
  Mnm="Ricker_Discrete.stan"
  Ndisc=length(unique(X))-1# # of gamma effects to estimate. Don't need to estimate X[iyr]=0
  if(Ndisc==1) Mnm="Ricker_Discrete_1level.stan"
  model_data<-list(Nyrs=Nyrs,SP=SP,mu_obslgRS=mu_obslgRS,sd_obslgRS=sd_obslgRS,X=X,Ndisc=Ndisc)
  
  ini_alpha=as.double(reg$coefficients[1]);ini_beta=as.double(reg$coefficients[2])
  if(ini_beta>0) ini_beta=-0.0001
  
  ini_gamma=rep(0,Ndisc)
  
  
  inits1=list(alpha=ini_alpha,beta=ini_beta,gamma=ini_gamma,sd_pro=1)
  inits2=list(alpha=ini_alpha*1.5,beta=ini_beta*0.75,gamma=ini_gamma,sd_pro=1.25)
  inits3=list(alpha=ini_alpha*0.75,beta=ini_beta*1.25,gamma=ini_gamma,sd_pro=0.75)
  
} else {
  Mnm="Ricker.stan"#assumess error around lgR/S is normally distributed
  ini_alpha=as.double(reg$coefficients[1]);ini_beta=as.double(reg$coefficients[2])
  if(ini_beta>0) ini_beta=-0.0001
  
  inits1=list(alpha=ini_alpha,beta=ini_beta,gamma=ini_gamma,sd_pro=1)
  inits2=list(alpha=ini_alpha*1.5,beta=ini_beta*0.75,gamma=ini_gamma*0.8,sd_pro=1.25)
  inits3=list(alpha=ini_alpha*0.75,beta=ini_beta*1.25,gamma=ini_gamma*1.2,sd_pro=0.75)
}

inits<-list(inits1,inits2,inits3)

#Call model and save fit object
fit=stan(file=Mnm,data=model_data,init=inits,chains=nchains,iter=niter,include=T)
Snm=paste0(OutDir,"Fit_",DoSite,"_",Am,"_",CoVarNm,".Rdata"); save(fit,file=Snm) 
