library("rstan")
#library("mnormt")
library("mixtools")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

logit<-function(x){
  log(x/(1-x))
}

inv_logit<-function(x){
  exp(x)/(1+exp(x))
}

RunMainstem=T
#OutDir="Output/NoCovEffect/"
#par(mfrow=c(2,2),mai=c(0.5,0.3,0.7,0.3), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
par(mfrow=c(2,1),mai=c(1,1,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
#par(mfrow=c(1,1),mai=c(1.2,0.8,0.2,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.3,cex.lab=1.2,font.lab=2)
for(itype in 1:1){
  OutDir=switch(itype,"Output/NoCovEffect/","Output/WithCovEffect/")
  for(itr in 1:2){#1:4
    if(RunMainstem==F){
        doTrib=switch(itr,"ubc","lcc")
    } else {
      doTrib=switch(itr,"knights landing","tisdale")
    }
    source ("GetData.R")
    fnm=paste0(OutDir,doTrib,".Rdata")
    load(file=fnm)
    
    dphi=as.data.frame(fit, pars = c("phi"))
    dlam=as.data.frame(fit, pars = c("lambda"))
    dmu=as.data.frame(fit, pars = c("mu_phi","mu_lambda","muRT","sigmaRT","rho"))
    
    phi=colMeans(dphi);lambda=colMeans(dlam) #year-specific meands
    #For MVN plot
    muRT=vector(length=2);sigmaRT=muRT
    icol=which(names(dmu)=="muRT[1]");muRT[1]=mean(dmu[,icol])
    icol=which(names(dmu)=="muRT[2]");muRT[2]=mean(dmu[,icol])
    icol=which(names(dmu)=="sigmaRT[1]");sigmaRT[1]=mean(dmu[,icol])
    icol=which(names(dmu)=="sigmaRT[2]");sigmaRT[2]=mean(dmu[,icol])
    icol=which(names(dmu)=="rho");rho=mean(dmu[,icol])
    vcv=matrix(nrow=2,ncol=2)
    vcv[1,1] = mean(sigmaRT[1]^2)
    vcv[2,2] = mean(sigmaRT[2]^2)
    vcv[1,2] = mean(rho*sigmaRT[1]*sigmaRT[2])
    vcv[2,1] = vcv[1,2]
    
    #points for confidence limit of MVN (alpha). ellipse returns x,y  positions at the specified
    #probability (alpha = 0.8 =80% CI). As they are in logit (x) and log (y) space, convert them
    
    xx=ellipse(muRT, vcv, alpha = 0.2, npoints = 250, newplot = F,draw = F)#x and y positions at specified alpha
    conx=inv_logit(xx[,1]);cony=exp(xx[,2])
    ymin=min(c(lambda,cony));ymax=max(c(lambda,cony))
    #xmin=min(c(phi,conx));xmax=max(c(phi,conx))*1.05
    
    
    
    if(RunMainstem==F){
      xmin=0.3;xmax=0.5
      ymin=10
      ymax=switch(itr,250,150)
    } else {
      xmin=switch(itr,0.375,0.30)
      xmax=switch(itr,0.475,0.55)
      ymin=switch(itr,20,10)
      ymax=switch(itr,70,120)
    }
    #if(itype==1){
    #xmin=switch(itr,0.35,0.3); xmax=switch(itr,0.45,0.5)
    #} else {
    #  xmin=switch(itr,0.3,0.35);xmax=switch(itr,0.5,0.55)
    #}
    xticks=seq(from=xmin,to=xmax,by=0.05)
    
    plot(phi,lambda,type='p',pch=19,cex=1.2,axes=F,bty='n',xlim=c(xmin,xmax),ylim=c(ymin,ymax),main=doTrib,xlab="",ylab="")
    axis(2);axis(1,at=xticks,labels=xticks)
    for(j in 1:length(xticks)){
      irec=min(which(theta>xticks[j]))
      axis(3,at=xticks[j],labels=calwk[irec],cex=0.8,pos=ymax*0.8)
      #axis(1,at=xticks[j],labels=paste0(calwk[irec],"        "),cex=1,outer=F,las=3)
    }
   
    #across-year means
    points(x=mean(dmu$mu_phi),mean(dmu$mu_lambda),type='p',pch=3,cex=1.8,col="red")#plot means of mvn
    
    #MVN countours
    lines(x=conx,y=cony,type='l',lwd=2,col='red') #80%
    
    xx=ellipse(muRT, vcv, alpha = 0.5, npoints = 250, newplot = F,draw = F)#x and y positions at specified alpha
    conx=inv_logit(xx[,1]);cony=exp(xx[,2])
    lines(x=conx,y=cony,type='l',lty=2,lwd=2,col='red')
    
    
    xyoff=1#1.01#label year-specific values
    for(iyr in 1:Nyrs){
      yr= paste0("'",substr(year[iyr],start=3,stop=4))
      text(x=phi[iyr]*xyoff,y=lambda[iyr]*xyoff,labels=yr,cex=0.9,pos=3)
    }
    
  }
}
mtext("Phi (median run date)",side = 1,line = 1, outer = T,cex=1.3,font=2)
mtext("Lambda (run timing steepness)",side=2, line=1,las=3,outer=T,cex=1.3,font=2)
