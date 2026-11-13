#Get data for model

fn_dat=paste0(OutDir,"SR_",DoSite,"_",Am,".dat")
d=read.table(file=fn_dat,header=T)

if(UseCovar==F){
  Nyrs=dim(d)[1]
  sryr_ind=1:Nyrs
  Nyrdrop=0
  
} else { #using a covariate so might drop some SR years if there is no covariate value
  
  if(UseSQcovar==F){#use all covariate data years in fitting
    fn_cov=paste0(OutDir,"Cov_",DoSite,"_",Am,".dat")
    
  } else { #only use years where all covariates are available (for model comparison)
    fn_cov=paste0(OutDir,"CovSQ_",DoSite,"_",Am,".dat")
  }
  
  dCov0=read.table(file=fn_cov,header=T)
  if(is.na(CoVarNm)==F){#when running a covariate model
    dCov=subset(dCov0,Variable==CoVarNm) #dCov will be used to set X below
    
  } else { #when running no covariate model but wanting to use same set of years as covariate models
    unvars=unique(dCov0$Variable)# a trick to get # of years for CovarNm=NA
    dCov=subset(dCov0,Variable==unvars[1]) #dCov will not be used to set X. but needed beolow to set Nyrs
  }
  
  #Identify the indices for stock-recruit years with covariate data
  sryr_ind=which(is.na(match(d$Year,dCov$Year))==F)
  Nyrs=length(sryr_ind)
  Nyrdrop=dim(d)[1]-Nyrs
}

SP=vector(length=Nyrs);R=SP;Byr=SP
mu_obslgR=SP;sd_obslgR=SP
mu_obslgRS=SP;sd_obslgRS=SP
X=rep(0,Nyrs);X0=X

for(iyr in 1:Nyrs){
  
  irow=sryr_ind[iyr]#the row in stock-recruit table for this iyr which has both SR and covariate data.
  
  Byr[iyr]=d$Year[irow]
  SP[iyr]=d$Stock[irow]
  if(SP[iyr]==0){
    SP[iyr]=1
    print(c(iyr, " spawners = 0. Assuming spawners=1"))
  }
  #************************************
  #Outmigrant abundance (recruitment)
  R[iyr]=d$muJuv[irow]*1000#convert from '000s to standard units (1=1 outmigrant)
  #*******************************
  mu_obslgR[iyr]=log(R[iyr])
  sd_obslgR[iyr]=sqrt(log(d$cvJuv[irow]^2+1))  #Calculated sd on log scale based on cv on linear scale

  #Simulate variation in lgR and divide by spawner stock to get mean and sd of obslgR/S
  simR=exp(rnorm(n=1000,mean=mu_obslgR[iyr],sd=sd_obslgR[iyr]))
  obslgRS=log(simR/SP[iyr])
  mu_obslgRS[iyr]=mean(obslgRS)
  sd_obslgRS[iyr]=sd(obslgRS)
  
  if(UseCovar==T){
    
    if(is.na(CoVarNm)==F){
      dCov1=subset(dCov,Year==d$Year[irow]) #Grab covariate data associated with current year in SR table
      X0[iyr]=dCov1$Value#raw
      X[iyr]=(X0[iyr]-mean(dCov$Value))/sd(dCov$Value)#standardized
    
      } else {
        X0[iyr]=0;X[iyr=0]
    }
  }
 
}
#if(Nyrdrop>0) print(c("Yrs dropped due to missing covariate values ",Nyrdrop))

