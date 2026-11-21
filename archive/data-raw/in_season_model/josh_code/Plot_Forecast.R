library("rstan")
library("scales")
rm(list = ls())

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

RunMainstem=T
DoForecast=T

#lq=0.16;uq=0.84#lower and upper quantiles for 1 sd below and above mean (central 66% of posterior)
lq=0.1;uq=0.9
#lq=0.025;uq=0.975
#OutDir="Output/NoCovEffect/"
OutDir="Output/lg1NoCovEffect/"

#Setup data to demonstrate error in forecast. A future year with mean and sd of abundance estimates
#through Dec, Jan, Feb, and Mar

#Units of N from spring run files are already in units of '000s.
xlims=c(0.2,0.6) #range of x-axis

#For_ewk_lab=c("Dec 30", "Jan 29","Feb 26","Apr 01")
For_ewk=c(17,22,26,30) #sequential week starting (Sep-02 = ewk 1)
Nfor=length(For_ewk)


par(mfrow=c(2,2),mai=c(0.6,0.5,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
if(DoForecast==F)par(mfcol=c(2,2),mai=c(0.6,0.5,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
#par(mfrow=c(2,1),mai=c(1,1,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)

for(itr in 1:4){#1:4
  doTrib=switch(itr,"ubc","lcc","knights landing","tisdale")
  RunMainstem=switch(itr,F,F,T,T)
  
  source ("GetData.R")
  fnm=paste0(OutDir,doTrib,".Rdata")
  load(file=fnm)
  
  #get the mean and sd of abundance across forecast weeks for current trib
  For_Nx_mu=vector(length=Nfor);For_Nx_sd=For_Nx_mu
  for(ifor in 1:Nfor){
    irow=For_ewk[ifor]
    Nxmu=as.double(Nx_mu[irow,]);Nxsd=as.double(Nx_sd[irow,]) #convert matrix to vector (one col for each yr)
    irecs=which(Nxmu!=0)#identify years where sample ended by this forecast week
    For_Nx_mu[ifor]=mean(Nxmu[irecs]) #remove these years from the stats
    For_Nx_sd[ifor]=mean(Nxsd[irecs])
  }
  
  df=as.data.frame(fit, pars = c("For_cp"))
  niter=dim(df)[1]
  hyp_stats=matrix(nrow=Nwks,ncol=3)
  for(iwk in 1:Nwks){hyp_stats[iwk,1:3]=quantile(df[,iwk],probs=c(lq,0.5,uq),na.rm=T)}
  
  #Plot forecasted run timing and abundance
  plot(theta,hyp_stats[,2],xlim=xlims,type='l',lwd=2,bty='n',,main=doTrib,xlab="",ylab="")#xlab="Proportion of Run Year" ylab="Run Timing Forecast"
  #polygon(x = c(theta, theta), y = c(hyp_stats[,1],hyp_stats[,3]), col = alpha("grey", 0.5), border = NA)
  polygon(x = c(rev(theta), theta), y = c(rev(hyp_stats[,1]),hyp_stats[,3]), col = alpha("grey", 0.5), border = NA)
  abline(v=For_ewk/Nwks,lty=2)
  par(srt=90)
  for(ifor in 1:Nfor){
    text(x=For_ewk[ifor]/Nwks,y=0.1,labels=calwk[For_ewk[ifor]],cex=0.9,pos=3)
  }
  par(srt=0)
  

  if(DoForecast==T){
    #Make a total abundance forecast based on Nx for each of the forecast weeks
    stats=matrix(nrow=Nfor,ncol=3)
    cvfor=vector(length=Nfor)
    for(ifor in 1:Nfor){
      p=df[,For_ewk[ifor]]
      
      #Uncertainty in cummulative abundance on forecast week
      pNx=rnorm(n=niter,mean=For_Nx_mu[ifor],sd=For_Nx_sd[ifor])
      
      #Calculate annual abundance by expanding vummulative abundance on forecast week
      pNtot=exp(pNx)/p        # by forcasted proportion of run that has passed
      cvfor[ifor]=sd(pNtot)/mean(pNtot)
      stats[ifor,1:3]=quantile(pNtot,probs=c(lq,0.5,uq),na.rm=T)
    }
    #ymin=min(stats)
    ymin=0
    ymax=max(stats[,2])
    plot(1:Nfor,stats[,2],pch=19,cex=1.2,bty='n',xlim=c(0.75,Nfor+0.25),ylim=c(ymin,ymax),axes=F,xlab="",ylab="",main=doTrib)#xlab="ForeCast Date",ylab="Forecasted Abundance"
    axis(2);axis(1,at=1:Nfor,labels=calwk[For_ewk])
    for(ifor in 1:Nfor){
      arrows(x0=ifor,x1=ifor,y0=stats[ifor,1],y1=stats[ifor,3],angle=90,code=3,length=0.025)
      text(x=ifor,y=ymax*0.95,labels=round(exp(For_Nx_mu[ifor]),digits=1),cex=0.9)
      #text(x=ifor,y=ymax*0.85,labels=round(For_Nx_sd[ifor],digits=2),cex=0.9)
      #text(x=ifor,y=ymax*0.75,labels=round(cvfor[ifor],digits=2),cex=0.9)
        
    }  
    
    if(itr==2 | itr==4){
      
      mtext("Proportion of year                                                                          ",side = 1,line = 1, outer = T,cex=1.3,font=2)
      mtext("Proportion of annual outmigrant abundance (0-1)",side=2, line=1,las=3,outer=T,cex=1.3,font=2)
      mtext("                                                                                 Week used to make forecast",side = 1,line = 1, outer = T,cex=1.3,font=2)
      mtext("Forecasted annual abundance ('000s)",side=2, line=-27,las=3,outer=T,cex=1.3,font=2)
      
    }
  }
}#next itr

mtext("Proportion of year                                                                          ",side = 1,line = 1, outer = T,cex=1.3,font=2)
mtext("Proportion of annual outmigrant abundance (0-1)",side=2, line=1,las=3,outer=T,cex=1.3,font=2)
#mtext("                                                                                 Week used to make forecast",side = 1,line = 1, outer = T,cex=1.3,font=2)
#mtext("Forecasted annual abundance ('000s)",side=2, line=-27,las=3,outer=T,cex=1.3,font=2)

