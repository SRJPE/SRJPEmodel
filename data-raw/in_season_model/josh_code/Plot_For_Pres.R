
library("rstan")
library("scales")


graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

RunMainstem=T
PlotObsError=F
OutDir="Output/NoCovEffect/"
#OutDir="Output/WithCovEffect/"
#OutDir="Output/lg1NoCovEffect/"
#OutDir="Output/lg1WithCovEffect/"

xlims=c(0.2,0.8) #better control of x axis and labelling
#xlims=c(0,1)
xvals=seq(from=xlims[1],to=xlims[2],by=0.2)

#doTrib="ubc"
doTrib="knights landing"
doYr=16

source ("GetData.R")
fnm=paste0(OutDir,doTrib,".Rdata")
load(file=fnm)

dp=as.data.frame(fit, pars = c("phi","lambda","muRT","bCov","cp"))

df=as.data.frame(fit, pars = c("For_cp"))
hyp_med=vector(length=Nwks)
for(iwk in 1:Nwks)hyp_med[iwk]=median(df[,iwk],na.rm=T)

#Uncertainty in weekly cummulative abundance
iyr=doYr
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

par(mfrow=c(1,1),mai=c(1.2,0.8,0.2,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.3,cex.lab=1.2,font.lab=2)
xmain=paste0(year[iyr]," (",round(mean(Nx_cv[,iyr]),digits=2),")")#Ntot_cv[iyr]
plot(pwk[1:Nestwks[iyr],iyr],exp(Nx_mu[1:Nestwks[iyr],iyr]),xlim=xlims,ylim=c(0,ymax),main=paste0(doTrib," ", year[doYr]),type='p',pch=19,cex=1.2,bty='n',axes=F,xlab="",ylab="")
axis(1,at=xvals,labels=xvals);axis(2)
for(j in 1:length(xvals)){
  irec=min(which(theta+0.001>=xvals[j]))
  axis(1,at=xvals[j],labels=paste0(calwk[irec],"        "),cex=1,outer=F,las=3)
}
legend("bottomright",legend=c("BT-SPAS-X","beta","beta+error","hyper"),pch=c(19,NA,NA,NA),lty=c(NA,1,1,1),lwd=c(NA,1,1,2),col=c("black","red","black","blue"),bty='n')
mtext("Proportion of year / date",side = 1,line = 1, outer = T,cex=1.3,font=2)
mtext("Cummulative annual outmigrant abundance",side=2, line=1,las=3,outer=T,cex=1.3,font=2)


#plot uncertainty in cummulative weekly abundances
if(PlotObsError==T){
  for(iwk in 1:Nestwks[iyr]) arrows(x0=pwk[iwk,iyr],x1=pwk[iwk,iyr],y0=statsNx[iwk,1],y1=statsNx[iwk,2],col="grey",angle=90,code=3,length=0.025)
}


#plot beta without deviates. Will include any modelled covariate effects as they are embedded in phi or lambda
icol=which(names(dp)==paste0("phi[",iyr,"]"));phi=dp[,icol]#phi=mean(dp[,icol])
icol=which(names(dp)==paste0("lambda[",iyr,"]"));lambda=dp[,icol]#lambda=mean(dp[,icol])
a=lambda*phi;b=lambda*(1-phi)

xstats=matrix(nrow=Nwks,ncol=3)
for(iwk in 1:Nwks){
  y_beta <- exp(Ntot_mu[iyr])*pbeta(theta[iwk], shape1 = a, shape2 = b) 
  xstats[iwk,]=as.double(quantile(y_beta,probs=c(0.025,0.5,0.975)))
  #lines(x=c(theta[iwk],theta[iwk]),y=c(xstats[iwk,1],xstats[iwk,3]),col="red")
}
lines(theta,xstats[,2],lty=1,lwd=1,col="red")

stats=matrix(nrow=Nwks,ncol=3)
for(iwk in 1:Nwks){
  icol=which(names(dp)==paste0("cp[",iwk,",",iyr,"]"))
  cp=dp[,icol]*exp(Ntot_mu[iyr])
  stats[iwk,]=as.double(quantile(cp,probs=c(0.025,0.5,0.975),na.rm=T))
}
lines(theta,stats[,2],lwd=2)
polygon(x = c(theta, theta), y = c(stats[,1],stats[,3]), col = alpha("grey", 0.5), border = NA)

#The beta run timing from hyper to compare with annual beta's.  
lines(theta,hyp_med*exp(Ntot_mu[iyr]),lty=1,lwd=2,col="blue")
