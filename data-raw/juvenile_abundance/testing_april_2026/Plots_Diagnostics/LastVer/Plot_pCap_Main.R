library(rstan)
library(dplyr)
library(scales)
library(brms)

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

inv_logit<-function(x){
  exp(x)/(1+exp(x))
}

lb=0.025;ub=0.975	#0.975

d0a=SRJPEdata::years_to_include_rst_data
d0=subset(d0a,site=="knights landing"| site=="tisdale")
uniqsite=unique(d0$site);Nsites=length(uniqsite)

par(mfrow=c(2,1),mai=c(1,1,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)

for(isite in 1:Nsites){
  d2=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]

  #observed trap efficiences for each trial
  pCap_inputs <- prepare_pCap_inputs(mainstem =T,mainstem_site=uniqsite[isite])
  Releases=pCap_inputs$inputs$data$Releases
  Recaptures=pCap_inputs$inputs$data$Recaptures
  obs_pCap=Recaptures/Releases
  flow=pCap_inputs$inputs$data$mr_flow

  load(paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_mainstem_skew_re_",DoSite,".Rdata"))
  dp=as.data.frame(pcap,pars=c("logit_pCap","b0_pCap","b_flow","pro_sd_P","alpha"))

  #Predicted trap efficiences for each trial
  Nobs=length(obs_pCap)
  pstats=matrix(nrow=Nobs,ncol=3)
  Vnm="logit_pCap"
  for(i in 1:Nobs){
    icol=which(names(dp)==paste0(Vnm,"[",i,"]"))
    pstats[i,]=quantile(inv_logit(dp[,icol]),probs=c(lb,0.5,ub));pstats[i,2]=mean(dp[,icol])
    pstats[i,2]=mean(inv_logit(dp[,icol]))
  }
  print(c(DoSite,round(cor(obs_pCap,pstats[,2]),digits=2)))

  #flow-trap efficiency relationship stats
  Nq=50;qmin=min(flow)*1.1;qmax=max(flow)*1.1;qseq=seq(qmin,qmax,length.out=Nq)
  Vnm="b0_pCap";icol=which(names(dp)==Vnm);b0_pCap=dp[,icol]
  Vnm="b_flow";icol=which(names(dp)==Vnm); b_flow=dp[,icol]
  Vnm="pro_sd_P";icol=which(names(dp)==Vnm); pro_sd_P=dp[,icol]
  Vnm="alpha";icol=which(names(dp)==Vnm); alpha=dp[,icol]
  Vnm="yr_sd_P";dp=as.data.frame(pcap,pars=Vnm);yr_sd_P=dp[,1]

  stats=matrix(nrow=Nq,ncol=3);statswe=stats
  Nsims=length(pro_sd_P)
  pdev=rskew_normal(n=Nsims, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)

  for(i in 1:Nq){
    pred=inv_logit(b0_pCap+b_flow*qseq[i])
    stats[i,]=quantile(pred,probs=c(lb,0.5,ub))#;stats[i,]=mean(pred)

    pred2=inv_logit(b0_pCap+b_flow*qseq[i] + pdev)
    statswe[i,]=quantile(pred2,probs=c(lb,0.5,ub))#;statswe[i,2]=mean(pred2)
  }

  #Plot it all
  xmin=min(c(flow,qmin));xmax=max(c(flow,qmax))
  ymin=min(c(obs_pCap,stats,pstats));ymax=max(c(obs_pCap,stats,pstats))
  plot(flow,obs_pCap,xlim=c(xmin,xmax),ylim=c(ymin,ymax),main=DoSite,pch=21,cex=1.5,col='blue',bty='l',xlab="",ylab="")

  points(flow,pstats[,2],pch=19,cex=1,col="black")
  for(i in 1:Nobs)lines(x=c(flow[i],flow[i]),y=c(pstats[i,1],pstats[i,3]),lwd=1)

  polygon(x = c(rev(qseq), qseq), y = c(rev(stats[,1]), stats[,3]), col = alpha("grey", 0.6), border = NA)
  polygon(x = c(rev(qseq), qseq), y = c(rev(statswe[,1]), statswe[,3]), col = alpha("grey", 0.3), border = NA)
  lines(qseq,stats[,2],lwd=2)

}

mtext("Standardized Discharge",side = 1,line = -1, outer = T,cex=1.1,font=2)
mtext("Capture Probability",side=2, las=3,outer=T,cex=1.1,font=2,line=-1)

legend("topright",legend=c("observed","predicted","predicted relationsip","95% CI (no process error)","95% CI (with process error)"),
       pch=c(21,19,NA,15,15),pt.cex=c(1.5,1,NA,1.5,1.5),col=c("blue","black","black",alpha("grey", 0.6),alpha("grey", 0.3)),
       lwd=c(NA,NA,1,NA,NA),bty='n')



#histograms showing distribution of observed and predicted, where latter uses mean and random year and process error
tp=0.25;xcol=c(rgb(1,0,0,tp),rgb(0,0,1,tp))
for(isite in 1:Nsites){
  d2=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]

  #observed trap efficiences for each trial
  pCap_inputs <- prepare_pCap_inputs(mainstem =T,mainstem_site=uniqsite[isite])
  Releases=pCap_inputs$inputs$data$Releases
  Recaptures=pCap_inputs$inputs$data$Recaptures
  obs_pCap=Recaptures/Releases


  load(paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_mainstem_skew_re_",DoSite,".Rdata"))
  dp=as.data.frame(pcap)

  Vnm="b0_pCap";icol=which(names(dp)==Vnm);b0_pCap=dp[,icol]
  Vnm="pro_sd_P";icol=which(names(dp)==Vnm); pro_sd_P=dp[,icol]
  Vnm="alpha";icol=which(names(dp)==Vnm); alpha=dp[,icol]
  Vnm="yr_sd_P";dp=as.data.frame(pcap,pars=Vnm);yr_sd_P=dp[,1]

  Nsims=dim(dp)[1]
  yrdev=rnorm(n=Nsims,mean=0,sd=yr_sd_P)
  pdev=rskew_normal(n=Nsims, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)
  genP_wyr=inv_logit(b0_pCap + yrdev + pdev)

  xmax=0.05 #limit plotted distributions to this value so can see fit in meat of distribution
  irecs1=which(obs_pCap<=xmax);irecs2=which(genP_wyr<=xmax)
  xlims=c(0,max(c(obs_pCap[irecs1],genP_wyr[irecs2])))
  xbreaks=seq(from=0,to=max(xlims),length.out=40)

  usevar_mean=mean(genP_wyr) #don't use truncated data to calculate mean
  usevar=genP_wyr[irecs2]

  pd1=hist(obs_pCap[irecs1],plot=F,breaks=xbreaks)$density
  pd2=hist(usevar,plot=F,breaks=xbreaks)$density

  ylims=c(0,max(c(pd1,pd2)))
  hist(obs_pCap[irecs1], xlab="",main=DoSite,ylab="",xlim=xlims,col=xcol[1],breaks=xbreaks,freq=F,axes=F,ylim=ylims)
  abline(v=mean(obs_pCap),col="red",lty=1,lwd=2)
  axis(1)
  hist(usevar,add=2,col=xcol[2],breaks=xbreaks,freq=F,ylim=ylims)
  abline(v=usevar_mean,col="blue",lty=2,lwd=2)

}
legend("topright",legend=c("Observed","Predicted"),col=xcol,pch=c(15,15),pt.cex=c(2,2),bty='n')
mtext("Capture Probability",side = 1,line = -1, outer = T,cex=1.1,font=2)

