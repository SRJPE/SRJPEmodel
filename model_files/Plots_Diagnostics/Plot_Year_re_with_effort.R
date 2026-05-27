rm(list = ls())
library(tidyverse)
library(SRJPEdata)
library(SRJPEmodel) # and any others
library(rstan)
#library(viridisLite)
library("brms")


graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,height=14,width=18);par(las=1)

logit<-function(x){log(x/(1-x))}
inv_logit<-function(x){exp(x)/(1+exp(x))}

#use_pdev=1 #includes process error. Black points show in a year with MR data for weeks without
use_pdev=0 #exclude process error. Easier to see variation among years.
use_effort=0

#IsMain=T;MainSite="tisdale"
#IsMain=T;MainSite="knights landing"
IsMain=F;Mainsite=NA

if(IsMain==F){#The site and year associated with each site-year random effect
  pCap_inputs <- prepare_pCap_inputs(model_type = "all_sites")

  Ntribs=pCap_inputs$inputs$data$Ntribs
  xlabels=paste(strtrim(pCap_inputs$site_year_fit$site,3),pCap_inputs$site_year_fit$run_year,sep="-")
  load("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_all_sites.Rdata")
  parlist=c("b0_pCap","pro_sd_P","yr_sd_P")#,"yr_re")#,"alpha")
  dp=as.data.frame(pcap,pars=c("b0_pCap","yr_re"," yr_sd_P","pro_sd_P"))

} else {
  pCap_inputs <- prepare_pCap_inputs(model_type = "one_site", skew = T, site_selection = MainSite)
  xlabels=pCap_inputs$years_fit

  fn=paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/pCap_one_site_skew_re_",MainSite,".Rdata")
  load(fn)
  parlist=c("pro_sd_P","yr_sd_P")#,"yr_re","alpha")
  dp=as.data.frame(pcap,pars=c("b0_pCap","yr_re"," yr_sd_P","pro_sd_P","alpha"))
}
Nsims=dim(dp)[1]
#print(summary(pcap,pars=parlist)$summary)

obs_pCap=pCap_inputs$inputs$data$Recaptures/pCap_inputs$inputs$data$Releases
ind_trib=pCap_inputs$inputs$data$ind_trib

lb=0.025;ub=0.975

#####Plot all random year effects on one plot ##############
par(mfcol=c(1,1),mai=c(0.85,0.85,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
Vnm="yr_re"

#Nyr_re=pCap_inputs$inputs$data$Nyr_re
Nyr_re=pCap_inputs$inputs$data$Nyr_re

Vnm="yr_re"
yr_re=matrix(nrow=Nyr_re,ncol=3)
for(i in 1:Nyr_re){
  vn=paste0(Vnm,"[",i,"]");icol=which(names(dp)==vn)
  yr_re[i,]=as.double(quantile(dp[,icol],prob=c(lb,0.5,ub)))
}


ylims=c(min(yr_re),max(yr_re))
plot(1:Nyr_re, yr_re[,2],type='p',pch=19,ylim=ylims,xlab="",ylab="Random Year Effect",bty='l',axes=F),legend=c("bottomright",legend=pCap_inputs$sites_fit)
for(i in 1:Nyr_re)  arrows(x0=i,x1=i,y0=yr_re[i,1],y1=yr_re[i,3],angle=90,code=3,length=0.025)

abline(h=0,lty=2)
axis(2)
axis(1,at=1:Nyr_re,labels=xlabels,las=3,cex.axis=0.8)
abline(v=1:Nyr_re,lty=3,col="gray")

if(IsMain==F){
  #Color-code points by site
  xlb=strtrim(xlabels,width=3);xlb_uniq=unique(xlb);nxlb=length(xlb_uniq)
  xcols=c(rep(c("black","red","yellow","purple","pink","green"),2),"black")#xcols=rainbow(n=nxlb)  #xcols=viridis(n=nxlb,option="D")

  pt_col=vector(length=Nyr_re)
  for(i in 1:nxlb){
    irecs=which(xlb==xlb_uniq[i])
    pt_col[irecs]=xcols[i]
  }
  for(i in 1:Nyr_re) points(i,yr_re[i,2],pch=19,col=pt_col[i])
}


### Plot site-year pCaps #############

if(IsMain==T){
  par(mfcol=c(1,1),mai=c(0.85,0.85,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
} else {
  par(mfrow=c(4,4),mai=c(0.5,0.6,0.1,0.1), omi=c(0.1,0.1,0.1,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
}


pCap_stats=matrix(data=0,nrow=Nyr_re,ncol=3)
for(i in 1:Nyr_re){

  vnm=paste0("yr_re[",i,"]");icol=which(names(dp)==vnm);yr_re=dp[,icol]

  if(IsMain==F){

    j=which(pCap_inputs$sites_fit==pCap_inputs$site_year_fit$site[i])
    vnm=paste0("b0_pCap[",j,"]");icol=which(names(dp)==vnm);b0_pCap=dp[,icol]
    Vnm=paste0("pro_sd_P[",j,"]");icol=which(names(dp)==Vnm);pro_sd_P=dp[,icol]
    Vnm=paste0("yr_sd_P[",j,"]");icol=which(names(dp)==Vnm);yr_sd_P=dp[,icol]
    pdev=rnorm(n=Nsims,mean=0,sd=pro_sd_P)

  } else {

    vnm="b0_pCap";icol=which(names(dp)==vnm);b0_pCap=dp[,icol]
    Vnm="pro_sd_P";icol=which(names(dp)==Vnm);pro_sd_P=dp[,icol]
    Vnm="alpha";icol=which(names(dp)==Vnm);alpha=dp[,icol]
    Vnm="yr_sd_P";icol=which(names(dp)==Vnm);yr_sd_P=dp[,icol]
    pdev=rskew_normal(n=Nsims, mu = 0, sigma = pro_sd_P, alpha = alpha, xi = NULL, omega = NULL)
  }
  yrdev=rnorm(n=Nsims,mean=0,sd=yr_sd_P)

  irecs=which(pCap_inputs$inputs$data$ind_yr==i)#get records of effort associated with this year
  effort=pCap_inputs$inputs$data$effort[irecs]
  p1=inv_logit(b0_pCap + yr_re + pdev*use_pdev + use_effort*log(mean(effort)))

  pCap_stats[i,]=as.double(quantile(p1,prob=c(lb,0.5,ub)));pCap_stats[i,2]=mean(p1)
  #print(c(i,mean(effort)))#to show how effort varies across years 1:Nyr_re
}



if(IsMain==T){

  #This is mean to show what error would be in a year without MR data
  #p1=inv_logit(b0_pCap + yrdev + pdev*use_pdev + use_effort*log(mean(mu_yr_effort)))
  effort=pCap_inputs$inputs$data$effort
  forP = inv_logit(b0_pCap + yrdev + pdev*use_pdev + use_effort*log(mean(effort)))
  forP=as.double(quantile(p1,probs=c(lb,0.5,ub)));forP[2]=mean(forP); forPmedian=median(forP)

  ylims=c(min(c(pCap_stats,forP)),max(c(pCap_stats,forP)))
  plot(1:Nyr_re,pCap_stats[,2],pch=19,cex=1.2,bty='l',xlim=c(1,(Nyr_re+1)),ylim=ylims,main=MainSite,axes=F,xlab="",ylab="")
  if(use_pdev==1) abline(h=mean(obs_pCap),lty=3,col="blue",lwd=3)

  abline(h=forP[2],lty=3,col="red")
  #abline(h=forPmedian,lty=2,col="purple")

  for(j in 1:Nyr_re){
    arrows(x0=j,x1=j,y0=pCap_stats[j,1],y1=pCap_stats[j,3],angle=90,code=3,length=0.025,lwd=1)
  }
  axis(2)
  axis(1,at=1:(Nyr_re+1),labels=c(pCap_inputs$years_fit,"no_et"),las=3)

  points(Nyr_re+1,forP[2],pch=19,cex=1.2,col="red")
  arrows(x0=Nyr_re+1,x1=Nyr_re+1,y0=forP[1],y1=forP[3],angle=90,code=3,length=0.025,col="red")

} else {


  for(i in 1:Ntribs){

    #pCap in a year without an estimated random effect
    vnm=paste0("b0_pCap[",i,"]");icol=which(names(dp)==vnm);b0_pCap=dp[,icol]
    vnm=paste0("pro_sd_P[",i,"]");icol=which(names(dp)==vnm);pro_sd_P=dp[,icol]
    vnm=paste0("yr_sd_P[",i,"]");icol=which(names(dp)==vnm);yr_sd_P=dp[,icol]

    yrdev=rnorm(n=Nsims,mean=0,sd=yr_sd_P)
    pdev=rnorm(n=Nsims,mean=0,sd=pro_sd_P)

    irecs=which(pCap_inputs$inputs$data$sd_yr_ind==i)
    nyrs=length(irecs)

    irecs2=which(pCap_inputs$inputs$data$ind_trib==i)
    effort=pCap_inputs$inputs$data$effort[irecs2]
    forP = inv_logit(b0_pCap + yrdev + pdev*use_pdev + use_effort*log(mean(effort)))
    forP=as.double(quantile(forP,probs=c(lb,0.5,ub)));forP[2]=mean(forP); forPmedian=median(forP)

    print(c(i,length(effort),mean(effort)))
    ylims=c(min(c(pCap_stats[irecs,],forP)),max(c(pCap_stats[irecs,],forP)))
    xmain=paste0(pCap_inputs$site_year_fit$site[irecs[1]]," n= ",length(irecs2)," yr_sd= ", round(mean(yr_sd_P),digits=2)," pro_sd= ", round(mean(pro_sd_P),digits=2),")")


    plot(1:nyrs,pCap_stats[irecs,2],pch=19,cex=1.2,bty='l',xlim=c(1,(nyrs+1)),ylim=ylims,main=xmain,axes=F,xlab="",ylab="")

    for(j in 1:nyrs)  arrows(x0=j,x1=j,y0=pCap_stats[irecs[j],1],y1=pCap_stats[irecs[j],3],angle=90,code=3,length=0.025)
    axis(2)
    axis(1,at=1:(nyrs+1),labels=c(pCap_inputs$site_year_fit$run_year[irecs],"no_et"),las=3)

    irecs2=which(ind_trib==i)
    if(use_pdev==1) abline(h=mean(obs_pCap[irecs2]),lty=3,col="blue",lwd=3)
    abline(h=forP[2],lty=3,col="red")
    #abline(h=forPmedian,lty=2,col="purple")

    points(nyrs+1,forP[2],pch=19,col="red")
    arrows(x0=nyrs+1,x1=nyrs+1,y0=forP[1],y1=forP[3],angle=90,code=3,length=0.025,col="red")
  }
}
mtext("Run Year",side = 1,line = -1, outer = T,cex=1.1,font=2)
if(use_effort==1){
  mtext("Effort Adjusted Capture Probability",side=2, las=3,outer=T,cex=1.1,font=2,line=-1)
} else{
  mtext("Capture Probability",side=2, las=3,outer=T,cex=1.1,font=2,line=-1)
}

