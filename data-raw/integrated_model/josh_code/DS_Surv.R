#1) Get the median size of spring run outmigrants by model week
#2) Get the multi-year proportion of spring run outmigrants leaving by model week from inseason model
#3) Get survival from RST to Delta for each model week given fish size for that week (1)
#4) Calculate a weighted average survival across run, with 2) providing the weights.

set.seed(123)

par(mfrow=c(2,3),xaxs="i",mai=c(.8,.7,.5,.5), omi=c(0.25,0.25,0.5,0.25),cex.axis=0.9)
PlotType=1 #model week vs. pRun and PLAD-based median size
#PlotType=2 #histogram of comp_post vs simulated as is done in JPE.stan 

#SurvCov=1
#Load CJS survival model fit object. 
if(SurvCov==0){
  fn_surv="fit_CovIndWY_FL.Rdata"#null model so only one element (iwy=1) for second dimension in SurvForecastSz and TribSurvForecastSz
  iwy=1# only one index in SurvForecastSz and TribSurvForecastSz for WY dimentsion
} else {
  fn_surv="fit_CovIndWY3_FL.Rdata"
  iwy=SurvCov #Here SurCov represents the element # to use for forecast year (e.g., iwy=1 for critical, iwy=2 for D/BN/AN, iwy=3 for W)
}


load(fn_surv);
#dsv=as.data.frame(fit,pars=c("SurvForecast"));print(c(mean(dsv[,1]),sd(dsv[,1])))
dsv=as.data.frame(fit,pars=c("SurvForecastSz","TribSurvForecastSz"))#SurvForecastSz[ix,1] #Only one WY class (all types) in current fit object (24 size classes)
Nsims=dim(dsv)[1]

SurvTypes=3 #Types of survival predictions from CJS model (mainstem, butte, feather)
TribStype=c(1,1,1,1,2,3)#ubc, ucc, mill, deer = 1, butte =2, and eventually feather =3

#Read in CJS forecasted survival posteriors
Nsz=25
Xsz  <- seq(from=10,to=130,length.out=Nsz)
SurvSz_stats=array(data=NA,dim=c(SurvTypes,Nsz,3))#First 3 is for mainstem (1 = ubc, ucc, mill, deer), butte (2, okie dam), feather (3, not inclueded yet)
SurvPost=array(dim=c(Nsims,Nsz,SurvTypes))

VnmM="SurvForecastSz";VnmT="TribSurvForecastSz"
for(isz in 1:Nsz){
  #mainstem
  Vnm=paste0(VnmM,"[",isz,",",iwy,"]")#1 is the WY class. For this version of model, WY was not included as a covarate so everything has index =11
  icol=which(names(dsv)==Vnm)
  jtrib=1
  SurvPost[1:Nsims,isz,jtrib]=dsv[,icol]
  SurvSz_stats[jtrib,isz,1:3]=as.double(quantile(dsv[,icol],probs=CI));SurvSz_stats[jtrib,isz,2]=mean(dsv[,icol])
  
  for(j in 1:(SurvTypes-1)){#butte, feather
    Vnm=paste0(VnmT,"[",isz,",",iwy,",",j,"]")
    icol=which(names(dsv)==Vnm)
    jtrib=j+1
    SurvPost[1:Nsims,isz,jtrib]=dsv[,icol]
    SurvSz_stats[jtrib,isz,1:3]=as.double(quantile(dsv[,icol],probs=CI));SurvSz_stats[jtrib,isz,2]=mean(dsv[,icol])
  }
  
}

Fl_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year
Surv_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year from weighted posterior sample
Surv_mu_check=vector(length=Ntribs)#pRun_weighted average survival calculated intuitive way for checking

#mean and sd of logit tranformed posterior samples of survival and the pRun-weighted mean forklength for JPE.stan model
DS_surv_mu=vector(length=Ntribs);DS_surv_sd=DS_surv_mu 

for(itrib in 1:Ntribs){
  
  # 1) Get across-year mean size of spring run by week from PLAD  #########
  #First step is to load the PLAD sizes for each trib and week
  fn_sz=paste0(PLADpath,"/SRforklength_",DoTrib[itrib],".csv")
  dsz=read.csv(file=fn_sz,header=T)
  Jwk=as.integer(rownames(dsz))
  Nwks=length(Jwk)# # of wks where sizes are available from PLAT
  
  # 2) Get proportion of run passing trap for every week in year (Sep-Aug) ########
  if(itrib==1){
    pRun=matrix(nrow=Ntribs,ncol=Nwks)
    Sz=pRun;  Surv=pRun
    
    labwk=character(length=Nwks)#for plotting
    for(iwk in 1:Nwks)labwk[iwk]=calwk[which(mwk==Jwk[iwk])]
  }
  Sz[itrib,1:Nwks]=dsz$X0.5# put size in a matrix for later plotting 
  
  #Load inseason outmigrant timing model to calculate proportion of run for each week
  fnm=paste0(InseasDir,DoTrib[itrib],"_null.Rdata")
  load(file=fnm)
  dcp=as.data.frame(fit, pars = c("For_cp"))
  
  iwk=1;imwk=which(mwk==Jwk[iwk])#identify the inseaon model week for the Julian week designation in PLAD
  pRun[itrib,iwk]=mean(dcp[,imwk],na.rm=T)
  for(iwk in 2:Nwks){
    imwk=which(mwk==Jwk[iwk])
    pRun[itrib,iwk]=mean(dcp[,imwk],na.rm=T)-sum(pRun[itrib,1:(iwk-1)])
  }
  pRun[itrib,]=pRun[itrib,]/sum(pRun[itrib,])#Make sure it sums to one to avoid any discretization error

  
  #3) Get pRun-weighted size and survival rate for outmigrant run ################
  
  #Set the jtrib index needed for SurvForecastSz. TribSurvForecastSz, and SurvSz_stats
  if(DoTrib[itrib]=="ubc"| DoTrib[itrib]=="ucc"| DoTrib[itrib]=="mill creek"| DoTrib[itrib]=="deer creek"){
    jtrib=1
  } else if (DoTrib[itrib]=="okie dam"){
    jtrib=2
  } else {#Feather eventually
    jtrib=3
  }
  
  Surv_mu_check[itrib]=0;Fl_mu[itrib=0]
  #Set of samples to grab from posterior sample of CJS survival for size class associated with following weeks
  if(itrib==1) isamps=sample(1:Nsims,size=500,replace=F) #to keep comp_post to reasonable size, sample 100 posterior survival values from each size class (wk)
  for(iwk in 1:Nwks){
    
    #find the closest size in Xsz given median size in this week from plad (Sz) to get the isz index for Xsz used in CJS model
    irecs=which(Xsz<=Sz[itrib,iwk])
    if(length(irecs)==0){ 
      isz=1 #Sz is smaller than lowest value of Xsz so set to lowest index
    } else{
      isz=max(irecs) #For Sz>max(Xsz) this will set index to last element for Xsz
    }

    #Create a composite posterior of survival rates which accumulates posterior samples across all weeks/size classes.
    #Posterior sample each wk/size class is weighted based on pRun relative to lowest pRun aross weeks (where reps=1)
    preps=round(pRun[itrib,iwk]/min(pRun[itrib,])) # # of replicates for this iwk
    if(iwk==1){
      comp_post=rep(SurvPost[isamps,isz,jtrib],times=preps) #repeat the random sample of posterior for this wk preps times
    } else {
      comp_post=c(comp_post,rep(SurvPost[isamps,isz,jtrib],times=preps))#each set of repeated samples is added to the vectory
    }
            
    Fl_mu[itrib]=Fl_mu[itrib]+pRun[itrib,iwk]*Sz[itrib,iwk]#pRun-weighted mean size
    
    #Most obvious way to compute a pRun-weighted average survival.
    Surv_mu_check[itrib]=Surv_mu_check[itrib]+pRun[itrib,iwk]*SurvSz_stats[jtrib,isz,2]
  }
  Surv_mu[itrib]=mean(comp_post)#To compare against Surv_mu_check
  
  #mean and sd of survival rate in logit space for outmigrant run for JPE.stan model
  DS_surv_mu[itrib]=mean(logit(comp_post))
  DS_surv_sd[itrib]=sd(logit(comp_post))
  
  print(c(DoTrib[itrib],round(Surv_mu_check[itrib],digits=4),round(Surv_mu[itrib],digits=4)))
  
  #Simulated survival rates as done in JPE.stan (but with lower constraint of -6.5 ~ 0.15% survival)
  #sim_surv used to see if composite posterior (comp_post) is accurately modelled by DS_surv_mu and DS_surv_sd and JPE.stan approap
  #For plotting (if PlotType==2) or text output on different plot type (if PlotType==1)
  sim_surv=inv_logit(rnorm(n=5000,mean=DS_surv_mu[itrib],sd=DS_surv_sd[itrib])) 
  
  if(PlotType==1){
    #Plot proportion of run and size by week with weighted survival mean and other stats printed on each panel
    plot(1:Nwks,pRun[itrib,],type='l',bty='n',main=DoTrib[itrib],xlab="",ylab="",axes=F)
    axis(2);axis(1,at=1:Nwks,labels=labwk,las=3)
    mu1=mean(sim_surv); sd1=sd(sim_surv)
    #text(x=1,y=max(pRun[itrib,])*0.9, labels=paste0("Wgt'd mean FL ",round(Fl_mu[itrib],digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.85, labels=paste0("Wgt'd mean Surv (Surv_mu_check) ",round(Surv_mu_check[itrib],digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.8, labels=paste0("Mean Surv from comp_post (Surv_mu) ",round(Surv_mu[itrib],digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.75, labels=paste0("mean of tranformed DSsurv ",round(mu1,digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.7, labels=paste0("sd Surv from comp_post ",round(sd(comp_post),digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.65, labels=paste0("sd of sim surv's from DSsurv_mu and _sd ",round(sd1,digits=2)),pos=4)
   
    
    par(new=T)
    plot(1:Nwks,Sz[itrib,],type='l',bty='n',col="red",xlab="",ylab="",axes=F)
    axis(4)
    if(itrib==Ntribs){
      legend("bottom",legend=c("proportion","forklength"),lty=c(1,1),col=c("black","red"),bty='n')
      mtext("Date",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
      mtext("Proportion of Outmigrant Run Passing Trap",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)
      mtext("Median Forklength of Spring Run Outmigrants (mm)",side=4, las=3,line=-1,outer=T,cex=1.3,font=2,col="red")
    }
    
  } else if (PlotType==2){
    Pr=density(sim_surv)
    Po=hist(comp_post,plot=F)  
    ymax=max(c(max(Po$density),max(Pr$y)))
    
    hist(comp_post,ylim=c(0,ymax),xlab="",ylab="",main=DoTrib[itrib],freq=F,axes=F)
    axis(1);lines(density(sim_surv),lty=2)
    
    if(itrib==Ntribs){
      legend("topright",legend=c("pRun-weighted composite posterior","JPE.stan simulated posterior"),lty=c(NA,2),pch=c(22,NA),pt.bg=c("gray","NA"),bty='n')
      mtext("Downstream Survival Rate",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
      mtext("Frequency",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)
    }
    
  }
}#next itrib


#Plot survival stats from CJS model and overlay stats that will be used in JPE.stan based on pRun-weighted mean size. Should be close
par(mfrow=c(2,2),xaxs="i",mai=c(.8,.7,.5,.5), omi=c(0.25,0.25,0.5,0.25),cex.axis=0.9)
xcol=c("black","blue","red","green")#multiple tribs use the same mainstem survival relationship
for(jtrib in 1:3){
  xmain=switch(jtrib,"ubc/ucc/mill/deer","Butte","Feather")
  plot(Xsz,SurvSz_stats[jtrib,,2],type='l',bty='n',ylim=c(0,max(SurvSz_stats[jtrib,,3])),xlab="",ylab="",main=xmain)
  polygon(x = c(rev(Xsz), Xsz), y = c(rev(SurvSz_stats[jtrib,,1]), SurvSz_stats[jtrib,,3]), col = alpha("grey", 0.25), border = NA)
  
  #check to see if simulated survival going into stan model (based on DS_surv_mu and DS_surv_sd
  #is reproducing what is coming out of CJS model 
  if(jtrib<=2){ #exclude feather for now
    TribStype=c(1,1,1,1,2,3)
    itribs=which(TribStype==jtrib)
    
    for(j in 1:length(itribs)){
    
      y1=inv_logit(rnorm(n=5000,mean=DS_surv_mu[itribs[j]],sd=DS_surv_sd[itribs[j]]))
      dstats=as.double(quantile(y1,probs=CI));dstats[2]=mean(y1)
      
      #For x-axis position, find the Xsz index closes to FL_mu, the pRun-weighted mean size 
      irecs=which(Xsz<=Fl_mu[itribs[j]])
      if(length(irecs)==0){isz_mu=1} else{isz_mu=max(irecs)}
      #points(Xsz[isz_mu],dstats[2],pch=19,cex=1.2,col=xcol[j])
      #lines(x=c(Xsz[isz_mu],Xsz[isz_mu]),y=c(dstats[1],dstats[3]),col=xcol[j])
    }
  }
  #if(jtrib==1) legend("topright",legend=c(DoTrib[1:4]),col=xcol,pch=rep(19,4),bty='n')
}
mtext("Forklength (mm)",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
mtext("Survival Rate to Delta",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)
