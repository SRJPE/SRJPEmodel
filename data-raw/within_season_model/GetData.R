
#ddate=read.table(file="../../RST/Data/Jwk_Dates.txt",header=T) #if you want to see calendar date and model index for each week
calwk=c('Sep-07','Sep-14','Sep-21','Sep-28','Oct-05','Oct-12','Oct-19','Oct-26','Nov-02','Nov-09',
        'Nov-16','Nov-23','Nov-30','Dec-07','Dec-14','Dec-21','Dec-28','Jan-04','Jan-11', 'Jan-18',
        'Jan-25','Feb-01','Feb-08','Feb-15','Feb-22','Mar-01','Mar-08','Mar-15','Mar-22','Mar-29',
        'Apr-05','Apr-12','Apr-19','Apr-26','May-03','May-10','May-17','May-24','May-31','Jun-07',
        'Jun-14','Jun-21','Jun-28','Jul-05','Jul-12','Jul-19','Jul-26','Aug-02','Aug-09','Aug-16',
        'Aug-23','Aug-30','Aug-31')

#Define period that will be modelled
fwk=36;lwk=35 #All year Sep-01 - Aug-31
mwk=c(fwk:53,1:lwk)#the calendar week associated with each modelled wk
Nwks=length(mwk) #the number of modelled weeks
theta=(1:Nwks)/Nwks;theta[Nwks]=0.999#the proportional week (proportion of total year through week iwk)

if(RunMainstem==F){
  RSTpath="../../RST/RunSize/StanVersion/Output_Trib_SR"
  RSTpathA="../../RST/RunSize/StanVersion/Output_Trib"
} else {
  RSTpath="../../RST/RunSize/StanVersion/Output_Main_SR"
  RSTpathA="../../RST/RunSize/StanVersion/Output_Main"
}

fn1=list.files(path=RSTpath,paste0("N_",doTrib))
Nyrs=length(fn1)
Nestwks=vector(length=Nyrs)#Number of weekly abundance estimates for each year
ewk=matrix(nrow=Nwks,ncol=Nyrs,data=0) #calender week for each weekly estimate
pwk=ewk #proportional week for each week-year combination. Varies by year based on when they sampled at RST

year=vector(length=Nyrs)
Ntot_mu=year;Ntot_sd=year;Ntot_cv=year
Nx_mu=matrix(nrow=Nwks,ncol=Nyrs,data=0)
Nx_sd=Nx_mu;Nx_cv=Nx_mu
CovX0=matrix(nrow=Nyrs,ncol=2)#annual covariates effecting phi [,1] and lambda [,2]
CovX=CovX0

#par(mfrow=c(4,4),mai=c(0.3,0.2,0.3,0.2), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
for(iyr in 1:Nyrs){
  year[iyr]=substr(fn1[iyr],start=nchar(fn1[iyr])-9,stop=nchar(fn1[iyr])-6)
  
  fn_fit=paste0(RSTpath,"/N_",doTrib,"_",year[iyr],".rdata")  
  load(fn_fit)
  
  #Ntot=abundance$sims.list$Ntot
  Ntot=abundance$Ntot  

  Ntot_cv[iyr]=sd(Ntot)/mean(Ntot) #for likelihood in log space
  Ntot_mu[iyr]=log(mean(Ntot))#original code
  #Ntot_mu[iyr]=mean(log(Ntot))#Alternative
  Ntot_sd[iyr]= sqrt(log(Ntot_cv[iyr]^2+1))
  
  fn_dat=paste0(RSTpathA,"/Inputs_",doTrib,"_",year[iyr],".rdata")  
  load(fn_dat)
  
  Jwk=abundance_inputs$weeks_fit
  Flow=abundance_inputs$catch_flow_raw
  Nstrata=abundance_inputs$inputs$data$Nstrata #u is in here (now called dd$u), as is Nstrata
  u=abundance_inputs$inputs$data$u
  
  TrimRun=F
  
  #strip out some of initial and last weeks with 0 catch
  if(TrimRun==T){
    wkbuff=3
    fuwk=min(which(u>0))-wkbuff;  if(fuwk<=0) fuwk=1
    luwk=max(which(u>0))+wkbuff;if(luwk>Nstrata) luwk=Nstrata
    Nestwks[iyr]=luwk-fuwk+1
  } else { #Use all weeks where abundance was estimated
    fuwk=1
    Nestwks[iyr]=Nstrata
  }

  for(iwk in 1:Nestwks[iyr]){
    ewk[iwk,iyr]=which(mwk==Jwk[iwk+fuwk-1])
    pwk[iwk,iyr]=theta[ewk[iwk,iyr]]
   
    icols=fuwk:(iwk+fuwk-1)#the
    
    #N=abundance$sims.list$N[,icols]
     N=abundance[,icols]#Assumes weekly estimates come first followed by Ntot
   
    if(iwk==1){
      cumN=N
    } else{
      cumN=rowSums(N) #abundance from start of sampling through week iwk
    }
    mu=mean(cumN)#the mean across posterior samples
    if(mu==0)  mu=.01 #can be 0 in cases with multiple weeks of zero catch at start of sampling
    Nx_mu[iwk,iyr]=log(mu)#total abundance through this week in log space
    sdx=sd(cumN)#sd across posterior samples of cummulative abundance
    if(sdx==0) sdx=0.01 #can also happen if multiple weeks of 0 catch at start of sampling
    Nx_cv[iwk,iyr]=sdx/mu#the cv in untranformed space
    Nx_sd[iwk,iyr]=sqrt(log(Nx_cv[iwk,iyr]^2+1))#convert from cv in untransformed space to sd in log space
    #Nx_sd[iwk,iyr]=sqrt(log(0.05^2+1)) #check to see if high cv is what is messing up fit
  }#iwk
  
  #Get annual flows statistic
  
  irecs=which(Jwk<=53 |(Jwk>=1 & Jwk<=5)) #Flows prior to Feb in yr 't'
  CovX0[iyr,1]=max(Flow[irecs])
  CovX0[iyr,2]=CovX0[iyr,1] #use same covariate for prediction of phi and lambda but this structure allows one to use different covariates

}#iyr

#Standardize annual covariate values
CovLabel=rep("Peak flow prior to February (cfs)",2)
for(j in 1:2){
  muX=mean(CovX0[,j]);sdX=sd(CovX0[,j])
  CovX[,j]=(CovX0[,j]-muX)/sdX
}

#mtext("Proportion of Year",side = 1,line = 1, outer = T,cex=1.3,font=2)
#mtext("Cummulative Abundance",side=2, line=1,las=3,outer=T,cex=1.3,font=2)


#dat=as.data.frame(cbind(year,Nestwks,Ntot_mu,Ntot_sd,t(Nx_mu),t(Nx_sd)))
#names(dat)=c("Year","NestWks","lgNtot_mu","lgNtot_sd",paste0("lgNx_mu",1:Nwks),paste0("lgNx_sd",1:Nwks))

dat=as.data.frame(cbind(year,Nestwks,Ntot_cv,t(Nx_cv)))
names(dat)=c("Year","NestWks","lgNtot_cv",paste0("lgNx_cv",1:Nwks))
fnout=paste0(OutDir,doTrib,".out")
write.table(file=fnout,x=dat,row.names=F,col.names=T)




