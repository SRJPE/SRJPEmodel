rm(list = ls())
library(scales)
library(readxl)
library("rstan")
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

logit<-function(x){log(x/(1-x))}
inv_logit<-function(x){exp(x)/(1+exp(x))}


RunModel=F#Run stan estimation model of annual abundance or read from saved fit objects
CI=c(0.1,0.5,0.9) #quantiles to compute stats on predictions

SRdir="c:/projects/BayDelta/SAC_JPE/stock_recruit/output/"
InseasDir="c:/projects/BayDelta/SAC_JPE/InSeason/BetaRT/output/rnd_error/"
PLADpath="c:/projects/BayDelta/SAC_JPE/RST/RunSize/StanVersion/PLAD_Size_Predictions"

calwk=c('Sep-07','Sep-14','Sep-21','Sep-28','Oct-05','Oct-12','Oct-19','Oct-26','Nov-02','Nov-09',
        'Nov-16','Nov-23','Nov-30','Dec-07','Dec-14','Dec-21','Dec-28','Jan-04','Jan-11', 'Jan-18',
        'Jan-25','Feb-01','Feb-08','Feb-15','Feb-22','Mar-01','Mar-08','Mar-15','Mar-22','Mar-29',
        'Apr-05','Apr-12','Apr-19','Apr-26','May-03','May-10','May-17','May-24','May-31','Jun-07',
        'Jun-14','Jun-21','Jun-28','Jul-05','Jul-12','Jul-19','Jul-26','Aug-02','Aug-09','Aug-16',
        'Aug-23','Aug-30','Aug-31')

mwk=c(36:53,1:35)#the calendar week associated with each inseason wk (Sep-Aug)
#bt_wk=c(45:53,1:22) #the weeks btspas abundance estimates were computed for (Nov-May) where we will have forklength predictions from PLAD

#dfor=read.csv(file="For_input.csv",header=T)#Read inputs to make forecast
dfor=read_excel("data-raw/integrated_model/josh_code/For_input.xlsx", sheet = "For_Wks",trim_ws=T,col_types=c("text"))
For_CalWk=dfor$For_CalWk
Nfor=length(For_CalWk)

dfor=read_excel("data-raw/integrated_model/josh_code/For_input.xlsx", sheet = "DS_Surv",trim_ws=T,col_types=c("numeric"))
SurvCov=dfor$SurvCov

dfor=read_excel("data-raw/integrated_model/josh_code/For_input.xlsx", sheet = "SR_InSeas",trim_ws=T,col_types=c("text","text","text","text",rep("numeric",times=2+Nfor*2)))
DoTrib=dfor$Trib;Am=dfor$Am;SRCovNm=dfor$SRCovNm;SRmodType=dfor$SRmodType
Sp=dfor$Sp  #currently using historical average
Xsr=dfor$Xsr_std#in standardized units, currently using mean across SQ years (thus std value =0)

Ntribs=length(DoTrib)
For_Nx_mu=matrix(nrow=Ntribs,ncol=Nfor);For_Nx_sd=For_Nx_mu

for(itrib in 1:Ntribs){
  ifor=1; For_Nx_mu[itrib,ifor]=dfor$For_Nx_mu_1[itrib];For_Nx_sd[itrib,ifor]=dfor$For_Nx_sd_1[itrib]
  ifor=2; For_Nx_mu[itrib,ifor]=dfor$For_Nx_mu_2[itrib];For_Nx_sd[itrib,ifor]=dfor$For_Nx_sd_2[itrib]
  ifor=3; For_Nx_mu[itrib,ifor]=dfor$For_Nx_mu_3[itrib];For_Nx_sd[itrib,ifor]=dfor$For_Nx_sd_3[itrib]
  ifor=4; For_Nx_mu[itrib,ifor]=dfor$For_Nx_mu_4[itrib];For_Nx_sd[itrib,ifor]=dfor$For_Nx_sd_4[itrib]
}

#These 3 scripts need to be run to produce parameters for JPE model (JPE.stan) and for plotting
source("SRfor.R")#get prior distribution of annual outmigrant abundance from SR model predicitons
source("InSeasfor.R") #get inseason-based prediction of annual outmigrant abundance (independent of SR priors)
source("DS_Surv.R")#get pRun-weighted downstream survival rate for each trib


### Run estimation model or read-in predictions from a previous run##
if(RunModel==T){
  source("Call_JPE.R")#stan model which estimates annual abundance using inseason data and priors from SR model
} else{
  fnm="Output/Fit.Rdata"#read in fit from previous run
  load(file=fnm)
}

pred_Ntot_stats=array(dim=c(Ntribs,Nfor,3))#Stats of posterior for forecast of annual outmigrant abundance
pred_DSsurv_stats=matrix(nrow=Ntribs,ncol=3)#Stats of posterior for forecast ofdownstream survival rate
JPE_trib_stats=array(dim=c(Ntribs,Nfor,3)) #Stats of posterior forForecast of JPE to Delta from each tributary

dp=as.data.frame(fit, pars = c("pred_Ntot","pred_Ntot_all","DS_surv","JPE_trib","JPE"))
Nsims=dim(dp)[1]
JPE_trib_post=array(dim=c(Nsims,Ntribs,Nfor))

for(itrib in 1:Ntribs){
  for(ifor in 1:Nfor){
    icol=which(names(dp)==paste0("pred_Ntot[",itrib,",",ifor,"]"));pred_Ntot=dp[,icol]
    pred_Ntot_stats[itrib,ifor, 1:3]=as.double(quantile(pred_Ntot,probs=CI))#;pred_Ntot_stats[itrib,ifor,2]=mean(pred_Ntot)

    icol=which(names(dp)==paste0("JPE_trib[",itrib,",",ifor,"]"));JPE_trib=dp[,icol]
    JPE_trib_stats[itrib,ifor, 1:3]=as.double(quantile(JPE_trib,probs=CI))
    JPE_trib_post[,itrib,ifor]=JPE_trib
  }

  icol=which(names(dp)==paste0("DS_surv[",itrib,"]"));DS_surv=dp[,icol]
  pred_DSsurv_stats[itrib, 1:3]=as.double(quantile(DS_surv,probs=CI))

}#next trib

pred_Ntot_all_stats=matrix(nrow=Nfor,ncol=4)
JPE_stats=matrix(nrow=Nfor,ncol=4)#Stats of posterior for JPE summed across tributaries
JPE_post=matrix(data=0,nrow=Nsims,ncol=Nfor)
for(ifor in 1:Nfor){
  icol=which(names(dp)==paste0("pred_Ntot_all[",ifor,"]"));pred_Ntot_all=dp[,icol]
  pred_Ntot_all_stats[ifor, 1:3]=as.double(quantile(pred_Ntot_all,probs=CI))
  pred_Ntot_all_stats[ifor,4]=sd(pred_Ntot_all)/mean(pred_Ntot_all)

  icol=which(names(dp)==paste0("JPE[",ifor,"]"));JPE=dp[,icol]
  JPE_stats[ifor, 1:3]=as.double(quantile(JPE,probs=CI))
  JPE_stats[ifor,4]=sd(JPE)/mean(JPE)
  JPE_post[,ifor]=JPE
}


###Some Plotting##
source("PlotFor.R") #plot SR and independent inseason forecasts for annual outmigrant abundance in comparison to model based ones (generated JPE.stan)
source("Plot_JPE.R") #details of predictions from JPE model (JPE.stan)


