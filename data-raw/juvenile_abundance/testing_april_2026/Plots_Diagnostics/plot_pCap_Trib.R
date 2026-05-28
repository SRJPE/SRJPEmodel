rm(list=ls(all=TRUE))

library(rstan)
library(scales)
library("brms")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)


inv_logit<-function(x){
  exp(x)/(1+exp(x))
}

lb=0.025;ub=0.975	#0.975

pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")
load("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_all_sites.Rdata")


options(max.print=5000)
parlist=c("b0_pCap","b_flow","pro_sd_P","yr_sd_P")#,"alpha")#,"yr_re")#,"alpha")
#print(summary(pcap,pars=parlist)$summary)#this is slow
#rhat = summary(pcap)$summary[,"Rhat"];ibad=which(rhat[1:length(rhat)-1]>1.05);print("# and vars with rhat>1.05");print(length(ibad));print(rhat[ibad])

#Get rid of inputs_site after running model again
Nmr=pCap_inputs$inputs$data$Nmr
Ntribs=pCap_inputs$inputs$data$Ntribs
ind_trib=pCap_inputs$inputs$data$ind_trib
mr_flow=pCap_inputs$inputs$data$mr_flow
Trib=pCap_inputs$sites_fit
Releases=pCap_inputs$inputs$data$Releases
Recaptures=pCap_inputs$inputs$data$Recaptures
pCap_obs=Recaptures/Releases

#plot predicted pCap for each trial vs 'observed' (r/R)
par(mfcol=c(1,1),mai=c(1,1,0.1,0.1), omi=c(0.1,0.1,0.3,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)

Vnm="logit_pCap";dp=as.data.frame(pcap,pars=Vnm)
pCap_pred=vector(length=length(pCap_obs))
for(i in 1:Nmr){
  vname=paste0(Vnm,"[",i,"]")
  icol=which(names(dp)==vname)
  pCap_pred[i]=mean(inv_logit(dp[,icol]))
}
xymin=min(pCap_obs,pCap_pred);xymax=max(pCap_obs,pCap_pred)
plot(pCap_obs,pCap_pred,type='p',pch=19,cex=0.7,xlim=c(xymin,xymax),ylim=c(xymin,xymax),xlab="Observed Capture Probability",ylab="Predicted Capture Probability",bty='l')
abline(coef=c(0,1),lty=2)
print("r2 for predicted vs observed pCaps across trials")
print(cor(pCap_obs,pCap_pred)^2)


#pCap distributions for each tributary in one panel and overlay hyper distribution
par(mfcol=c(1,1),mai=c(0.85,2.25,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
xmax=0.3
pstat=matrix(nrow=Ntribs,ncol=3)
dp=as.data.frame(pcap,pars=c("b0_pCap","pro_sd_P"))
Nsims=dim(dp)[1]

for(i in 1:Ntribs){
  vn=paste0("b0_pCap","[",i,"]");icol=which(names(dp)==vn);b0_pCap=dp[,icol]
  #vn=paste0("pro_sd_P","[",i,"]");icol=which(names(dp)==vn);pro_sd_P=dp[,icol]
  #vn=paste0("yr_sd_P","[",i,"]");icol=which(names(dp)==vn);yr_sd_P=dp[,icol]

  pC=inv_logit(b0_pCap)#+rnorm(n=Nsims,mean=0,sd=pro_sd_P))

  #use mean of fitted estimates for each trial
  #irecs=which(ind_trib==i)
  #pC=pCap_pred[irecs]
  pstat[i,]=as.double(quantile(pC,prob=c(lb,0.5,ub)));pstat[i,2]=mean(pC)
}

tnm=character(length=Ntribs);mup=vector(length=Ntribs)
mup=vector(length=Ntribs);mup2=mup
for(i in 1:Ntribs){
  irecs=which(ind_trib==i)
  nsamps=length(irecs)
  tnm[i]=paste0(Trib[i]," (",nsamps,")")
  recap=Recaptures[irecs]
  mup[i]=sum(Recaptures[irecs])/sum(Releases[irecs])
}

plot(pstat[,2],1:Ntribs,pch=19,bty='l',xlim=c(0,xmax),axes=F,ylab="",xlab="Mean Trap Efficiency")
axis(1);axis(2,at=1:Ntribs,labels=tnm,cex.axis=0.75)
for(i in 1:Ntribs){#plot CI's and point estimate for each trib
  arrows(x0=pstat[i,1],x1=pstat[i,3],y0=i,y1=i,angle=90,code=3,length=0.025,col="black")
  points(x=mup[i],y=i,pch=21,cex=1.5,col='red')
  #points(x=mup2[i],y=i,pch=21,cex=1.5,col='red')
}

Vnm="trib_mu_P";dp_mu=as.data.frame(pcap,pars=Vnm)#plot mean from hyper as vertical line1
hmu=inv_logit(mean(dp_mu[,1]));abline(v=hmu,lty=2)
ci=as.double(inv_logit(quantile(dp_mu[,1],probs=c(lb,ub))))
polygon(x=c(ci[1],ci[2],ci[2],ci[1]),y=c(1,1,Ntribs,Ntribs),col = alpha("darkgray", 0.3), border = NA)
#hmu2=sum(Recaptures)/sum(Releases);abline(v=hmu2,lty=2,col="blue")
#abline(v=mean(pstat[,2]),lty=2,col="red",lwd=3)

#hyper distribution
Vnm="trib_sd_P";dp_sd=as.data.frame(pcap,pars=Vnm)
pvec=seq(0,xmax,length.out=50);logit_pvec=log(pvec/(1-pvec))
prob=dnorm(x=logit_pvec,mean=mean(dp_mu[,1]),sd=mean(dp_sd[,1]))
par(new=T); plot(pvec,prob,lty=1,type='l',xlim=c(0,xmax),axes=F,xlab="",ylab="")
legend("topright",pch=c(19,21,NA),lty=c(NA,NA,1),col=c("black","red","black"),legend=c("Model Estimate","Point Estimate","Hyper-distribution"),bty='n')
########################################################################################


####Q-pCap relationships for each trib on separate plot along with data
par(mfrow=c(4,4),mai=c(0.5,0.6,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
xcol=rainbow(n=Ntribs)
Qs=seq(-2,6,by=0.1)
#Qs=seq(min(mr_flow),max(mr_flow),length.out=100)
Nq=length(Qs)
pstat=array(dim=c(Ntribs,Nq,3))
for(i in 1:Ntribs){
  vnm1=paste0("b0_pCap[",i,"]");dp1=as.data.frame(pcap,pars=vnm1);icol1=which(names(dp1)==vnm1)
  vnm2=paste0("b_flow[",i,"]");dp2=as.data.frame(pcap,pars=vnm2);icol2=which(names(dp2)==vnm2)
  for(j in 1:Nq){
    logit_pCap=dp1[,icol1]+dp2[,icol2]*Qs[j]#doing it for eac trial, which preserves correlation between intercept and slope
    tpCap=inv_logit(logit_pCap)
    pstat[i,j,1:3]=quantile(tpCap,probs=c(lb,0.5,ub))
    pstat[i,j,2]=mean(tpCap)
  }
}

for(i in 1:Ntribs){
  irecs=which(ind_trib==i)
  p=Recaptures[irecs]/Releases[irecs];Q=mr_flow[irecs]
  ylims=c(0,max(pstat[i,,],p))
  plot(Qs,pstat[i,,2],type='l',lwd=2,bty='l',ylim=ylims,col=xcol[i],xlab="",ylab="",main=Trib[i])#,xlim=c(min(Qs)*1.1,max(Qs)*1.2))
  polygon(x = c(rev(Qs), Qs), y = c(rev(pstat[i,,1]), pstat[i,,3]), col = alpha(xcol[i], 0.3), border = NA)
  points(Q,p,pch=19,cex=0.8,col=xcol[i])
  abline(v=0,lty=2)
  print(c(Trib[i],length(Q),round(mean(Q),digits=2)))
}
mtext("Standardized Discharge",side = 1,line = -1, outer = T,cex=1.1,font=2)
mtext("Capture Probability",side=2, las=3,outer=T,cex=1.1,font=2,line=-1)



#Compare distributions of pCap by trip including random year and process error effects
tp=0.25;xcol=c(rgb(1,0,0,tp),rgb(0,0,1,tp))
par(mfrow=c(4,4),mai=c(0.5,0.6,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
dp=as.data.frame(pcap,pars=c("b0_pCap","yr_sd_P","pro_sd_P"))#,"alpha"))
Nsims=dim(dp)[1]

for(i in 1:Ntribs){
  vn=paste0("b0_pCap","[",i,"]");icol=which(names(dp)==vn);b0_pCap=dp[,icol]
  vn=paste0("yr_sd_P[",i,"]");icol=which(names(dp)==vn);yr_sd_P=dp[,icol]
  vn=paste0("pro_sd_P[",i,"]");icol=which(names(dp)==vn);pro_sd_P=dp[,icol]
  #vn=paste0("alpha[",i,"]");icol=which(names(dp)==vn);alpha=dp[,icol]

  yrdev=rnorm(n=Nsims,mean=0,sd=yr_sd_P)
  #pdev=rskew_normal(n=Nsims, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)
  pdev=rnorm(n=Nsims,mean=0,sd=pro_sd_P)

  #pC=inv_logit(b0_pCap + yrdev + pdev)
  pC=inv_logit(b0_pCap +  pdev)

  irecs=which(ind_trib==i)
  pCap_o=pCap_obs[irecs]

  xmax=max(pCap_o)
  #if(xmax==0)  xmax=as.double(quantile(pC,prob=c(0.95)))

  irecs1=which(pCap_o<=xmax);irecs2=which(pC<=xmax)

  xlims=c(0,max(c(pCap_o[irecs1],pC[irecs2])))
  xbreaks=seq(from=0,to=max(xlims),length.out=40)

  pd1=hist(pCap_o[irecs1],plot=F,breaks=xbreaks)$density
  pd2=hist(pC[irecs2],plot=F,breaks=xbreaks)$density

  ylims=c(0,max(c(pd1,pd2)))

  hist(pCap_o[irecs1], xlab="",main=paste0(Trib[i]," (",length(irecs),")"),ylab="",xlim=xlims,col=xcol[1],breaks=xbreaks,freq=F,axes=F,ylim=ylims)
  abline(v=mean(pCap_o),col="red",lty=1,lwd=2)
  axis(1)
  hist(pC[irecs2],add=2,col=xcol[2],breaks=xbreaks,freq=F,ylim=ylims)
  abline(v=mean(pC),col="blue",lty=2,lwd=2)

  #hist(pC[irecs],main=Trib[i],xlab="",ylab="",axes=F);axis(1)
  #abline(v=median(pC),lty=1,lwd=2,col="red")
  #irecs=which(ind_trib==i)
  #abline(v=median(Recaptures[irecs]/Releases[irecs]),lwd=2,col="blue")
}
mtext("Capture Probability",side=1, las=1,outer=T,cex=1.1,font=2,line=-1)
legend("topright",legend=c("Observed","Predicted"),col=xcol,pch=c(15,15),pt.cex=c(2,2),bty='n')
