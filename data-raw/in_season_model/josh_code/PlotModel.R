
library("rstan")
library("scales")
library(SRJPEdata)
library(SRJPEmodel)

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)



PlotCovar=T
PlotObsError=F

#OutDir="Output/NoCovEffect/"
OutDir="Output/WithCovEffect/"
#OutDir="Output/lg1NoCovEffect/"
#OutDir="Output/lg1WithCovEffect/"

xlims=c(0.2,0.8) #better control of x axis and labelling
xvals=seq(from=xlims[1],to=xlims[2],by=0.2)

#par(mfrow=c(2,1),mai=c(1,1,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
par(mfrow=c(2,2),mai=c(1,1,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
for(itr in 1:4){#1:4
  #if(RunMainstem==F){
  #  doTrib=switch(itr,"ubc","lcc")
  #} else {
  #  doTrib=switch(itr,"knights landing","tisdale")
  #}
  
  doTrib=switch(itr,"ubc","lcc","knights landing","tisdale")
  RunMainstem=switch(itr,F,F,T,T)
  
  source ("GetData.R")
  
  fnm=paste0(OutDir,doTrib,".Rdata")
  load(file=fnm)
  
  dp=as.data.frame(fit, pars = c("phi","lambda","muRT","bCov","cp"))
  
  df=as.data.frame(fit, pars = c("For_cp"))
  hyp_med=vector(length=Nwks)
  for(iwk in 1:Nwks)hyp_med[iwk]=median(df[,iwk],na.rm=T)
  
  if(PlotCovar==F){
    if(itr==1){
      nrow=5;ncol=5
    } else if (itr==2 | itr==3){
      nrow=6;ncol=5
    } else {
      nrow=4;ncol=4
    }
    
    par(mfrow=c(nrow,ncol),mai=c(0.3,0.2,0.3,0.2), omi=c(1,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
    for(iyr in 1:Nyrs){
      #Uncertainty in weekly cummulative abundance
      statsNx=matrix(nrow=Nestwks[iyr],ncol=2)
      for(iwk in 1:Nestwks[iyr]){
        simNx=exp(rnorm(n=1000,mean=Nx_mu[iwk,iyr],sd=Nx_sd[iwk,iyr]))
        statsNx[iwk,]=as.double(quantile(simNx,probs=c(0.025,0.975)))
      }
      
      #Plot 'data' which are cummulative weekly abundances
      if(PlotObsError==F){
        ymax=max(exp(Nx_mu[,iyr]))*1.05
      } else {
        ymax=max(statsNx)*1.05
      }
      
      xmain=paste0(year[iyr]," (",round(mean(Nx_cv[,iyr]),digits=2),")")#Ntot_cv[iyr]
      plot(pwk[1:Nestwks[iyr],iyr],exp(Nx_mu[1:Nestwks[iyr],iyr]),xlim=xlims,ylim=c(0,ymax),type='p',pch=19,cex=1,bty='n',main=xmain,axes=F)
      axis(1,at=xvals,labels=xvals)
      if(iyr>(Nyrs-4)){#draw end of week data on bottom row of graphs
        for(j in 1:length(xvals)){
          irec=min(which(theta+0.001>=xvals[j]))
          axis(1,at=xvals[j],labels=calwk[irec],cex=0.8,outer=T,las=3)
        }
      }
      
      #plot uncertainty in cummulative weekly abundances
      if(PlotObsError==T){
        for(iwk in 1:Nestwks[iyr]) arrows(x0=pwk[iwk,iyr],x1=pwk[iwk,iyr],y0=statsNx[iwk,1],y1=statsNx[iwk,2],col="grey",angle=90,code=3,length=0.025)
        #polygon(x = c(pwk[1:Nestwks[iyr],iyr], pwk[1:Nestwks[iyr],iyr]), y = c(statsNx[,1],statsNx[,2]), col = alpha("pink", 0.5), border = NA)
      }
      
      stats=matrix(nrow=Nwks,ncol=3)
      for(iwk in 1:Nwks){
        icol=which(names(dp)==paste0("cp[",iwk,",",iyr,"]"))
        cp=dp[,icol]*exp(Ntot_mu[iyr])
        stats[iwk,]=as.double(quantile(cp,probs=c(0.025,0.5,0.975),na.rm=T))
      }
      lines(theta,stats[,2],lwd=2)
      polygon(x = c(theta, theta), y = c(stats[,1],stats[,3]), col = alpha("grey", 0.5), border = NA)
      
      #plot beta without deviates. Will include any modelled covariate effects as they are embedded in phi or lambda
      icol=which(names(dp)==paste0("phi[",iyr,"]"));phi=dp[,icol]#phi=mean(dp[,icol])
      icol=which(names(dp)==paste0("lambda[",iyr,"]"));lambda=dp[,icol]#lambda=mean(dp[,icol])
      a=lambda*phi;b=lambda*(1-phi)
      #y_beta <- exp(Ntot_mu[iyr])*pbeta(theta, shape1 = a, shape2 = b) 
      #lines(theta,y_beta,lty=1,lwd=2,col="red")
      xstats=matrix(nrow=Nwks,ncol=3)
      for(iwk in 1:Nwks){
        y_beta <- exp(Ntot_mu[iyr])*pbeta(theta[iwk], shape1 = a, shape2 = b) 
        xstats[iwk,]=as.double(quantile(y_beta,probs=c(0.025,0.5,0.975)))
        #lines(x=c(theta[iwk],theta[iwk]),y=c(xstats[iwk,1],xstats[iwk,3]),col="red")
      }
      lines(theta,xstats[,2],lty=1,lwd=1,col="red")
      #polygon(x = c(theta, theta), y = c(xstats[,1],xstats[,3]), col = alpha("red", 0.5), border = NA)
      
      #The beta run timing from hyper to compare with annual beta's.  
      lines(theta,hyp_med*exp(Ntot_mu[iyr]),lty=1,lwd=2,col="blue")
      
      if(iyr==2)   legend("bottomright",legend=c("BT-SPAS-X","beta+error","beta","hyper"),pch=c(19,NA,NA,NA),lty=c(NA,1,1,1),lwd=c(NA,1,1,2),col=c("black","black","red","blue"),bty='n')
    }
    mtext("Proportion of year / date",side = 1,line = 5, outer = T,cex=1.3,font=2)
    mtext("Proportion of annual outmigrant abundance (0-1)",side=2, line=1,las=3,outer=T,cex=1.3,font=2)
    mtext(doTrib,side = 3,line = 1, outer = T,cex=1.3,font=2)
  
    } else if (PlotCovar==T){
  
    #par(mfrow=c(1,1),mai=c(0.8,0.8,1,0.5), omi=c(0.8,0.8,0.8,0.8),cex.main=1.1,cex.lab=1.1,font.lab=2)
    #Plot covariate relationship
    icov=1
    xinc=seq(min(CovX[,icov]),max(CovX[,icov]),length.out=50)
    txinc=seq(min(CovX0[,icov]),max(CovX0[,icov]),length.out=50)
    icol=which(names(dp)==paste0("bCov[",icov,"]"));bCov=dp[,icol]
    icol=which(names(dp)==paste0("muRT[",icov,"]"));muRT=dp[,icol]
    stats=matrix(nrow=length(xinc),ncol=3)
    for (j in 1:length(xinc)){
      py=exp(muRT+bCov*xinc[j])/(1+exp(muRT+bCov*xinc[j]))
      stats[j,]=quantile(py,probs=c(0.025,0.5,0.975),na.rm=T)
    }
    ymin=0.30;ymax=0.55
    plot(txinc,stats[,2],type='l',lty=1,bty='n',axes=F,ylim=c(ymin,ymax),main=doTrib,xlab="",ylab="")#ylim=c(0.3,0.5)
    
    polygon(x = c(rev(txinc), txinc), y = c(rev(stats[,1]),stats[,3]), col = alpha("grey", 0.5), border = NA)
    axis(1)
    
    yticks=seq(from=ymin,to=ymax,by=0.05)
    axis(2,at=yticks,labels=yticks)
    for(j in 1:length(yticks)){
      irec=min(which(theta>yticks[j]))
      axis(4,at=yticks[j],labels=calwk[irec],cex=0.8,pos=min(txinc))
    }
    
  }
  
}
if(PlotCovar==T){
  mtext(CovLabel[icov],side = 1,line = -1, outer = T,cex=1.3,font=2)
  mtext("Median run date",side=2, line=-1,las=3,outer=T,cex=1.3,font=2)
}

  

