library("rstan")
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
nchains=3;niter=1000

MultiRun=T
Mnm="Ricker.stan"#assumess error around lgR/S is normally distributed

OutDir="Output/"


UseCovar=T#Covariates included in model except when unCovVarNm==NA
UseSQcovar=T#do you want to use all covariate data or just years where data are available for all covriates

for(itype in 1:4){
  DoSite=switch(itype,"ubc","ubc","lcc","lcc")
  Am=switch(itype,"rd","us","rd","us")
  
  fn_cov=paste0(OutDir,"CovSQ_",DoSite,"_",Am,".dat")
  dCov0=read.table(file=fn_cov,header=T)
  #NA to run model without a covariate but using years from square covariate table
  #Make sure it is at end so it uses the last Nyr which was based on a covariate (so uses reduced SR set)
  unCoVarNm=c(NA,unique(dCov0$Variable))
  Ncovars=length(unCoVarNm)
  
  for(irun in 1:Ncovars){
    CoVarNm=unCoVarNm[irun]
    print(c(irun,CoVarNm))
    source("Call_Model.R")
  }
  
}

