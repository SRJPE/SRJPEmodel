library(rstan)
library(SRJPEdata)

rm(list=ls(all=TRUE))

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

OutDir="Output/"
JuvDir="../RST/RunSize/StanVersion/Output_Main_SR/"

UseMill=F

DoSite="knights landing"#juvenile data
CovarStream="sacramento river"
Amethod="upstream_estimate";Am="us" #adult enumeration method in long and short form (the latter for file naming)

if(UseMill==T){
  DoStream=c("battle creek","clear creek","mill creek")#,"deer creek") #only 8 yrs of data for deer
} else {
  DoStream=c("battle creek","clear creek")
}

dA0=observed_adult_input #all adult data from this table in SRJPEdata
dA1=subset(dA0,is.na(match(stream,DoStream))==F & data_type==Amethod)

#Get list of years where escapement was measured in all DoStream s
dA2=subset(dA1,stream==DoStream[1]);yrlist1=unique(dA2$year)
dA2=subset(dA1,stream==DoStream[2]);yrlist2=unique(dA2$year)
if(UseMill==T){
  dA2=subset(dA1,stream==DoStream[3]);yrlist3=unique(dA2$year)
  Adult_yrs=Reduce(intersect, list(yrlist1, yrlist2, yrlist3))
  #dA2=subset(dA1,stream==DoStream[4]);yrlist4=unique(dA2$year);yrlist4=unique(dA2$year)
  #Adult_yrs=Reduce(intersect, list(yrlist1, yrlist2, yrlist3,yrlist4))
} else {
  Adult_yrs=Reduce(intersect, list(yrlist1, yrlist2))
}
dA=subset(dA1,is.na(match(year,Adult_yrs))==F)

#Files created from this script include full stock-recruit pairings given site and adult data type
fnout=paste0(OutDir,"SR_",DoSite,"_",Am,".dat") #SR data file
write(file=fnout,x="Year Stock muJuv cvJuv",ncolumns=1,append=F)
fn_cov=paste0(OutDir,"Cov_",DoSite,"_",Am,".dat") #will only include covariate years with a match for SR brood yer
fn_covSQ=paste0(OutDir,"CovSQ_",DoSite,"_",Am,".dat") #as above but no missing years for any of the covariates (square)

#Make SR table. Begin by getting list of years with outmigrant abundance estimates
file.ls=intersect(list.files(path=JuvDir,pattern=paste0("N_",DoSite)),list.files(path=JuvDir,pattern=".rdata"))
Nyrs=length(file.ls)

BroodYr=-99 #create a list of brood years with juvenile data. This first bunk record will be dropped at end of iyr loop

for(iyr in 1:Nyrs){
  JuvYr=as.integer(substr(x=file.ls[iyr],start=nchar(file.ls[iyr])-9,stop=nchar(file.ls[iyr])-6))
  
  #Check if there is an adult estimate for the Calendar year before the outmigrant run year
  #If so add a record to stock-recruit table
  Ayr=JuvYr-1
  irecs=which(dA$year==Ayr)
  print(c(Ayr,length(irecs)))
  if(length(irecs>0)){
    
    BroodYr=c(BroodYr,Ayr)
    
    # sum adult estimate across tribs for this year
    #Sp=sum(dA$passage_estimate[irecs]) #if using upstream_passage_estimates
    Sp=sum(dA$count[irecs])#if using observed_adult_input
   
    #Load total outmigrant abundance posterior
    fnjuv=paste0(JuvDir,file.ls[iyr]);load(fnjuv)
    dj=abundance$Ntot #In units of '000s. see PlAD_Integration.R in RST/RunSize/stan_version
    muJuv=mean(dj);cvJuv=sd(dj)/muJuv#get mean and CV
    
    write(file=fnout,x=c(Ayr,Sp,muJuv,cvJuv),ncolumns=4,append=T)
    
  } #has adult and juvenile data for this year
}#year loop
BroodYr=BroodYr[2:length(BroodYr)]#get rid of bunk record in first element


#Create covariate file which will only include BroodYr years (or fewer if no covariate for a brood year)
#given adult data of the type being selected. Missing covariate data for a brood year will be trapped in GetData.R
dCov0=stock_recruit_covariates
dCov=subset(dCov0,stream==CovarStream & site_group==DoSite & gage_number==11390500 & is.na(match(year,BroodYr))==F )
#dCov1=subset(dCov,is.na(value)==F & is.finite(value)==T)
#temporarily exclude temperature
dCov1=subset(dCov,is.na(value)==F & is.finite(value)==T)# & covariate_type=="flow")

#Create a new variable that combines lifestage and variable name
nrecs=dim(dCov1)[1];ls=character(length=nrecs)
irecs=which(dCov1$lifestage=="migration" | dCov1$lifestage=="spawning and incubation");ls[irecs]="mi"
irecs=which(dCov1$lifestage=="rearing");ls[irecs]="rr"
Variable=paste0(ls,"_",dCov1$covariate_structure)

dCov2a=as.data.frame(cbind(dCov1$year,Variable,dCov1$value));names(dCov2a)=c("Year","Variable","Value")

#add in 0 values to work with model types MeanlgRS.stan and Ricker_Discete.stan
nyr=length(BroodYr)
dCov2b= as.data.frame(cbind(BroodYr,rep("mu",nyr),rep(0,nyr)));names(dCov2b)=c("Year","Variable","Value")
dCov2c= as.data.frame(cbind(BroodYr,rep("null",nyr),rep(0,nyr)));names(dCov2c)=c("Year","Variable","Value")

#add in WY type index values
dWY0=read.csv(file="WY.csv",header=T)
irecs=which(is.na(match(dWY0$WY,BroodYr))==F)
dWY=dWY0[irecs,];nyr=length(irecs)

WY2ind=vector(length=nyr)
irecs=which(dWY$WYT=="C" | dWY$WYT=="D");WY2ind[irecs]=0
irecs=which(dWY$WYT=="BN");WY2ind[irecs]=1
irecs=which(dWY$WYT=="AN" | dWY$WYT=="W");WY2ind[irecs]=2

WY3ind=WY2ind
irecs=which(dWY$WYT=="C");WY3ind[irecs]=0
irecs=which(dWY$WYT=="D" |dWY$WYT=="BN");WY3ind[irecs]=1
irecs=which(dWY$WYT=="AN"| dWY$WYT=="W");WY3ind[irecs]=2

dCov2d=as.data.frame(cbind(BroodYr,rep("WY2",nyr),WY2ind));names(dCov2d)=c("Year","Variable","Value")
dCov2e=as.data.frame(cbind(BroodYr,rep("WY3",nyr),WY3ind));names(dCov2e)=c("Year","Variable","Value")
dCov2=rbind(dCov2b,dCov2c,dCov2d,dCov2e,dCov2a)

write.table(file=fn_cov,dCov2,row.names=F,quote=F)


#Finally, create a square covariate table where there are no missing years for all covariates (for model comparisons)
unyrs=unique(dCov2$Year)
unvars=unique(dCov2$Variable);nvars=length(unvars)
keepyr=""
for(iyr in 1:length(unyrs)){
  dCov3=subset(dCov2,Year==unyrs[iyr])
  nvars_for_yr=length(unique(dCov3$Variable))#dim(dCov3)[1]
  #print(c(unyrs[iyr],nvars_for_yr))
  if(nvars_for_yr==nvars) keepyr=c(keepyr,unyrs[iyr])
}
keepyr=keepyr[2:length(keepyr)]
irecs=which(is.na(match(dCov2$Year,keepyr))==F)
write.table(file=fn_covSQ,dCov2[irecs,],row.names=F,quote=F)


#Look at SR data
dSR=read.table(file=fnout,header=T)
plot(dSR$Stock,dSR$muJuv,bty='l',pch=19,cex=0,xlab="upstream passage",ylab="Outmigrant Abundance ('000s)")#,xlim=c(0,1500))
Nyrs=dim(dSR)[1]
for(iyr in 1:Nyrs){
  xlab=paste0("'",substr(dSR$Year[iyr],start=3,stop=4))
  text(x=dSR$Stock[iyr],y=dSR$muJuv[iyr],labels=xlab,cex=0.75)
}

