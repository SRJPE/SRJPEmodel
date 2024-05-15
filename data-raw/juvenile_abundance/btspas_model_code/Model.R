library("R2WinBUGS")
library(splines2)

EffortAdjust=T
MultiRun_Mode=T

if(MultiRun_Mode==F){
	#doTrib="deer creek_deer creek";doYr=2010
	doTrib="battle creek_ubc";doYr=2009
	#doTrib="clear creek_lcc";doYr=2005
	#doTrib="feather river_eye riffle";doYr=2020

	doWks=c(seq(45,53),seq(1,22)) #can't have a wider range than range used to build RST_input.txt (as specified in BuildData.R)
	#doWks=seq(1,19)

	fnprefix=paste0(doTrib,"_",doYr)
	#fnprefix=paste0("output/",doTrib,"_",doYr)
} else {
	fnprefix=paste0("data-raw/juvenile_abundance/btspas_model_code/OutSpecPriors/",doTrib,"_",doYr)
}

#Nmcmc=2500;Nburnin=500;Nthin=1;Nchains=3
#Nmcmc=5000;Nburnin=3000;Nthin=2;Nchains=3
Nmcmc=10000;Nburnin=5000;Nthin=5;Nchains=3
#Nmcmc=20000;Nburnin=10000;Nthin=10;Nchains=3

### Setup Data and initial values for pCap component of model ########
d0=read.table(file="data-raw/juvenile_abundance/RST_Input.txt",header=T)

#Tribs to use for pCap component of model. This is where sacramento river_knights landing is excluded because its efficiency is not exchangeable with those from CV tribs
alltribs=c("battle creek_ubc","clear creek_lcc","clear creek_ucc","feather river_eye riffle","feather river_herringer riffle","feather river_steep riffle","feather river_gateway riffle","feather river_sunset pumps")#,"sacramento river_knights landing"
Ntribs=length(alltribs)

#Trib index to use to predict efficiency for missing strata.
#Note length=0 if doTrib is not part of trib set that has MR data. In this case model that samples from trib hyper will be called
use_trib=which(alltribs==doTrib)

#Pick up data for current trib and run year, and exclude weeks that are not modelled.
dC=subset(d0,Trib==doTrib & RunYr==doYr & is.na(match(Week,doWks))==F)

#tribs to include in MR component of model (exchangeable). Used to be done using catch flow
dmr=subset(d0,is.na(match(Trib,alltribs))==F & is.na(Rel1)==F & is.na(Recap1)==F & is.na(mr_flow)==F)#catch_flow

# Unique MR experiments across tribs, years and weeks
mr_expts=unique(dmr[c("Trib","RunYr","Week")])
Nmr=dim(mr_expts)[1];allyrs=unique(mr_expts$RunYr)

Ntribs=length(alltribs);Nwks=length(doWks)
Releases=vector(length=Nmr);Recaptures=Releases;mr_flow=Releases
ind_trib=vector(length=Nmr);ind_yr=ind_trib;ind_wk=ind_trib

#Loop through MR set and assign variables and indices (what, trib, year and week it belongs to
#Note not all indices will be used in pCap model (e.g., ind_Wk won't be used if there is not a Week effect)
for(i in 1:Nmr){
	dmr1=subset(dmr,Trib==mr_expts$Trib[i] & RunYr==mr_expts$RunYr[i] & Week==mr_expts$Week[i])
	Releases[i]=dmr1$Rel1
	Recaptures[i]=dmr1$Recap1

	#create standardized flow for mark-recap experiments. catch_flow which is average for Jwk, mr_flow is average over recapture days (< 1 week)
	#use mr_flow in fitting (as for /pCap/pCap.tws), but standardized based on catch_flow, which will be used in prediction
	#for strata without MR data. Fitting based on mr_flow is more accurate and produces gentler b_flow.? slopes which provides more stability in abundance estimates
	#Otherwise you can get very low pCaps at high flows for some tris
	dmr1X=subset(dmr,Trib==mr_expts$Trib[i] & is.na(match(Week,doWks))==F)
	mr_flow[i]=(dmr1$mr_flow-mean(dmr1X$catch_flow,na.rm=T))/sd(dmr1X$catch_flow,na.rm=T)
	#mr_flow[i]=(dmr1$catch_flow-mean(dmr1X$catch_flow,na.rm=T))/sd(dmr1X$catch_flow,na.rm=T)

	ind_trib[i]=which(alltribs==mr_expts$Trib[i])
	ind_yr[i]=which(allyrs==mr_expts$RunYr[i])
	if(length(which(doWks==mr_expts$Week[i]))>0) ind_wk[i]=which(doWks==mr_expts$Week[i])
}

#Needed for PlotModel to show all MR data. Also good way to make sure indices are right.
MR=data.frame(cbind(mr_expts,ind_trib,ind_yr,ind_wk,Releases,Recaptures,mr_flow));names(MR)=c("Trib","RunYr","Wk","TribInd","YrInd","WkInd","Rel1","Recap1","mr_flow")
write.table(file="MRdata.txt",MR,col.names=T,row.names=F)

#######################################################################################################

#### Setup data and initial values for abundance component of model ##############
Nstrata=dim(dC)[1]
u0=dC$us	#may includes strata with no catch data which will be na

if(doTrib=="deer creek_deer creek" | doTrib=="mill creek_mill creek"){
	#currenly spring run not distingusihed, and there is also no effort for any records
	#for now use total chinook, and don't adjust as effort data not available
	u0=dC$ua
}

#Effort adjustment will be ratio of average effort for trib across all years and doWks period relative to weekly effort value
Effort=dC$Effort
if(length(which(is.na(Effort)==T))==Nstrata){#No effort data
	muEff=1	#Deer and mill have run years where there is no effort data, so assume constant over time
	Effort=1
} else {
	dCx=subset(d0,Trib==doTrib & is.na(match(Week,doWks))==F & is.na(Effort)==F)
	bad_recs=which(is.na(u0)==F & Effort==0) #In rare cases where strata was fished but no effort was recorded assume average effort for year
	if(length(bad_recs)==0){
		muEff=mean(dCx$Effort,na.rm=T)
	} else {
		muEff=mean(dCx$Effort[-bad_recs],na.rm=T)
		Effort[bad_recs]=muEff
	}
}
#Finally, adjust u based on ratio of mean effort to strata-specific effort
if(EffortAdjust==T){u=round(u0*muEff/Effort,digits=0) } else { u=u0 }

#Prior for upper limit on log abundance for any week. Note lgN estimated in units of thousands to avoid values in millions which occurs becasue of huge weekly catches in Butte
#Likely results in reduced precision problems and better mising for other systems with big numbers as well

#BT-SPAS constraint is 20 in log space but units are in 1's, not '000s.
#So for BT-SPAS-X convert 20, divide by 1000, then put back into log units.
#lgN_max=rep(log(exp(20)/1000),times=Nstrata)

#Default prior fornow
lgN_max=log((u/1000+1)/0.025)		#maxmim possible value for log N across strata
imiss=which(is.na(u)==T);lgN_max[imiss]=mean(lgN_max,na.rm=T) #for missing strata set to average max across strata

#lgN_max for special cases
dsp=read.csv(file="data-raw/juvenile_abundance/btspas_model_code/Special_Priors.csv",header=T)
dsp1=subset(dsp,Stream_Site==doTrib & RunYr==doYr)
istrata= which(is.na(match(doWks,dsp1$Jweek))==F)
lgN_max[istrata]=dsp1$lgN_max


#Calculate standardized flow for each weekly strata (can be different than mr_flow which only averages over days of recovery
#Note the mean and sd used for standardization has to be the same used for pCap model (above) for this trip
dX=subset(dC,Trib==doTrib & is.na(match(Week,doWks))==F)
cf=dC$catch_flow
irecs=which(is.na(dC$catch_flow)==T)
if(length(irecs)>0) cf[irecs]=mean(cf,na.rm=T) #if missing values for one or more strata set to mean (but no missing values as of Aug 02,2022.
catch_flow=(cf-mean(dX$catch_flow,na.rm=T))/sd(dX$catch_flow,na.rm=T)


#Identify the elements in 1:Nstrata (unmarked catch set) without and with pCap and corresponding flow data
#Alternate is to identify records without mr_flow but with catch_flow, and then sub catch_flow into mr_flow for these cases. Complicated and a bit manky.
Uind_woMR=which(is.na(dC$Rel1)==T | is.na(dC$mr_flow)==T);Nwomr=length(Uind_woMR);Jwks_womr=dC$Week[Uind_woMR]
Uind_wMR=which(is.na(dC$Rel1)==F & is.na(dC$mr_flow)==F);Nwmr=length(Uind_wMR);Jwks_wmr=dC$Week[Uind_wMR]
if(Nwomr==1) Uind_woMR=c(Uind_woMR,-99);if(Nwmr==1) Uind_wMR=c(Uind_wMR,-99)#so bugs doesn't bomb if only one strata with missing MR

#identify the elements in pCap from full MR dataset for U Strata being estimated
ind_pCap=which(MR$Trib==doTrib & MR$RunYr==doYr & is.na(match(MR$Wk,Jwks_wmr))==F)
Uwc_ind=which(is.na(u)==F);Nstrata_wc=length(Uwc_ind) #The elements of 1:Nstrata that have catch data (RST fished)

#Setup a data file for current run which lines up with output. Used for Plotmodel
trel=rep(NA,Nstrata);trel[Uind_wMR]=Releases[ind_pCap];trec=rep(NA,Nstrata);trec[Uind_wMR]=Recaptures[ind_pCap]
write.table(file=paste0(fnprefix,"_data.out"),cbind(dC$Week,u, trel,trec,round(Effort,digits=1),round(catch_flow,digits=1),round(dC$catch_flow,digits=1),round(exp(lgN_max),digits=1)),quote=F,row.names=F,col.names=c("Jwk","u","Releases","Recaptures","Effort","Std_Flow","Flow","N_max_000s"))

#### Setup B-spline basis matrix
k_int=4	#Rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)
Nknots=round(Nstrata/k_int,digits=0)
fkp=2;lkp=Nstrata-1 #Keep first and/or last knot positions away from tails if there are intervals with no sampling on the tails
kknots=seq(fkp,lkp,length.out=Nknots) #Define position of b-spline knots using even interval if no missing data
ZP <- bSpline(x=1:Nstrata,knots=kknots,deg=3,intercept=T)	 #bspline basis matrix. One row for each data point (1:Nstrata), and one column for each term in the cubic polynomial function (4) + number of knots
K=dim(ZP)[2]
#######################################################################################################

### Pass data and initial values to lists, define parameters to save, and run model ##################

#Select correct model given data for this run. Note slightly different data inputs across models
if(length(use_trib)==1){ #Some MR was done in this tributary

	if(Nwomr==0){#all strata have corresponding MR data
		ModName="estN_allMR"
		data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","ind_pCap","Nstrata_wc","Uwc_ind","mr_flow","lgN_max")

		#ModName="gamma_estN_allMR"
		#data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","ind_pCap","Nstrata_wc","Uwc_ind","mr_flow","lgN_max")

	} else {#some or all strata don't have MR data

		if(Nwmr>0){#some strata have MR data
			ModName="estN_missMR"		#MR data for some strata and tributary was included in the MR analysis

			data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","ind_pCap","Nwmr","Nwomr","Uind_wMR","Uind_woMR","use_trib","Nstrata_wc","Uwc_ind","mr_flow","catch_flow","lgN_max")

			#ModName="gamma_estN_missMR"		#MR data for some strata and tributary was included in the MR analysis
			#data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","ind_pCap","Nwmr","Nwomr","Uind_wMR","Uind_woMR","use_trib","Nstrata_wc","Uwc_ind","mr_flow","catch_flow","lgN_max")

		} else if (Nwmr==0){ #No strata have MR data
			ModName="estN_noMR"		#No MR data for any strata and trib was included in MR analysis
			data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","Nwomr","Uind_woMR","use_trib","Nstrata_wc","Uwc_ind","mr_flow","catch_flow","lgN_max")
		}
	}

} else if (length(use_trib)==0) { #No MR was done in this trib
	ModName="estN_noMR_notrib"	#No MR data for any strata and trib was not included in MR analysis
	data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","Nwomr","Uind_woMR","Nstrata_wc","Uwc_ind","mr_flow","catch_flow","lgN_max")#
}


#"beta_dev","tau.eta",
#parameters<-c("trib_mu.P","trib_sd.P","flow_mu.P","flow_sd.P","pro_sd.P","b0_pCap","b_flow","pCap_U","N","Ntot","sd.N","sd.Ne","sd.Ntot","pSplineVar")
parameters<-c("trib_mu.P","trib_sd.P","flow_mu.P","flow_sd.P","pro_sd.P","b0_pCap","b_flow","pCap_U","N","Ntot","sd.N","sd.Ne")

ini_b0_pCap=vector(length=Ntribs);for(itrib in 1:Ntribs){irows=which(ind_trib==itrib);ini_b0_pCap[itrib]=logit(sum(Recaptures[irows])/sum(Releases[irows]))}
ini_b0_pCap[4:8] <- c(NA, NA, NA, NA, NA)
#N is estimated in units of thousands in log space
ini_lgN=log(u/1000+2);irecs=which(is.na(ini_lgN)==T);ini_lgN[irecs]=log(2/1000)

inits1<-list(trib_mu.P=logit(sum(Recaptures)/sum(Releases)), b0_pCap=ini_b0_pCap, flow_mu.P=0, b_flow=rep(0,Ntribs),trib_tau.P=1, flow_tau.P=1,pro_tau.P=1, b_sp=rep(1,K),lg_N=ini_lgN)
inits2<-inits1;inits3<-inits1;inits<-list(inits1,inits2,inits3)


ModName2=here::here("data-raw", "juvenile_abundance", "btspas_model_code", paste0(ModName,".bug"))
post<-bugs(data, inits, parameters, ModName2,n.chains=Nchains, n.burnin=Nburnin,
           n.thin=Nthin,n.iter=Nmcmc, debug=F,codaPkg=F,DIC=T,clearWD=T,
           bugs.directory=here::here("data-raw", "WinBUGS14"))	#run bugs model

fnstats=paste0(fnprefix,"_post.out");write.table(post$sims.list,file=fnstats,col.names=T,row.names=T)
fnstats=paste0(fnprefix,"_sum.out");write.table(round(post$summary,digits=3),file=fnstats,col.names=T,row.names=T)
fn_dic=paste0(fnprefix,"_dic.out");write(file=fn_dic,c(post$pD,post$DIC),ncolumns=2)
fn_knots=paste0(fnprefix,"_knots.out");write(file=fn_knots,kknots,ncolumns=Nknots)

if(MultiRun_Mode==T){
	print(c(doTrib,doYr,ModName))
	Sys.sleep(0.01);flush.console()
}

