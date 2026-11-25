library("rstan")
library("loo") #https://mc-stan.org/loo/articles/loo2-with-rstan.html

rm(list=ls(all=TRUE))

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
nchains=3;niter=1000

OutDir="Output/"

#To run model without a covariate effect select CoVarNm=="null". This can be run with max SR years (UseSQCovar==F)
#or only SR years with values for all covariates (UseSCovar==T)

UseSQcovar=T#do you want to use all years with values for a given covariate,or just years where data are available for all covriates

for(itype in 5:5){
  DoSite=switch(itype,"ubc","ubc","lcc","lcc","mill creek","deer creek","okie dam","okie dam","knights landing")
  Am=switch(itype,"rd","us","rd","us","rd","hd","cs","hd","us")
  
  fn_cov=paste0(OutDir,"CovSQ_",DoSite,"_",Am,".dat")
  dCov0=read.table(file=fn_cov,header=T)
  
  unCoVarNm=c(unique(dCov0$Variable))
  Ncovars=length(unCoVarNm)

  n=vector(length=Ncovars);r2=n;sd_pro=n
  gamma_stats=matrix(nrow=Ncovars,ncol=4)
  looic_stats=matrix(nrow=Ncovars,ncol=4)#mean ic , se of ic, delta ic, and poverlap ic with best model
  
  for(irun in 1:Ncovars){#1:Ncovars
    
    CoVarNm=unCoVarNm[irun]
    
    print(c(irun,CoVarNm))
    source("Call_Model.R")              #fit model and save it to rdata file
    source("SumResults_FromMultiRun.R") #summary statistics and loo ic of model
  }
  
}

