d=read.csv(file=here::here("data-raw", "survival_model_STAN", "SacInp_withcov.csv"),header=T)

#Survival Reaches: 1=release-woodson, 2=woodson-butte, 3=butte-sac, 4=sac-delta
#Detection Locations: 1=release, 2=woodson, 3=butte, 4=sac, 5=delta
Nreaches=4
Ndetlocs=5
ReachKM=c(40, 88, 170, 110) #see Receiver_km sheet in Summary.xlsx

Nind=dim(d)[1]
Year=sort(unique(d$year));Nyrs=length(Year)

#Reorder release groups from order in file so they are ordered by year and release group in upstream-downstream direction
RelGp_raw=unique(d$StudyID); RelGp=RelGp_raw #Don't need to specify all RelGp[] values as in some cases RelGp_raw=RelGp (e.g., for 1)
RelGp[2]=RelGp_raw[11];RelGp[3]=RelGp_raw[12];RelGp[4]=RelGp_raw[2];RelGp[5]=RelGp_raw[3]
RelGp[6]=RelGp_raw[7];RelGp[7]=RelGp_raw[13];RelGp[10]=RelGp_raw[4];RelGp[11]=RelGp_raw[5];RelGp[12]=RelGp_raw[10];RelGp[13]=RelGp_raw[6]
Nrg=length(RelGp)

rch_covind=c(1,2,3,3)#index pointing to covariate effect for each reach (note Butte-Sac and Sac-Delta have same fixed effect index)

CH=matrix(nrow=Nind,ncol=Ndetlocs)
yrind=vector(length=Nind);rgind=yrind
WY2=yrind;WY3=yrind
firstCap=rep(0,Nind);lastCap=firstCap

for(i in 1:Nind){
  yrind[i]=which(Year==d$year[i])
  rgind[i]=which(RelGp==d$StudyID[i])

  #0 or 1 for dry (C,D,BN) and wet (W) water year types for dummmy variable effect
  if(d$year[i]==2013 | d$year[i]==2015 | d$year[i]==2020 | d$year[i]==2021){WY2[i]=0} else {WY2[i]=1}

  #3 water year type groupings (C, D-BN, W) to used as a fixed effect
  if(d$year[i]==2015 | d$year[i]==2021){
    WY3[i]=0
  } else if (d$year[i]==2013 | d$year[i]==2016 | d$year[i]==2018 | d$year[i]==2020){
    WY3[i]=1
  } else {
    WY3[i]=2
  }

  for(j in 1:Ndetlocs){
    CH[i,j]=as.integer(substr(d$ch[i],start=j,stop=j))
  }

  firstCap[i]=min(which(CH[i,1:Ndetlocs]==1))  #first station individual was detected at
  lastCap[i]=max(which(CH[i,1:Ndetlocs]==1))   #last station individual was detected at

}

#Data summary of releases and detections by release group and also do summary stats on individual covariates
FL=d$fish_length;WGT=d$fish_weight;CF=d$fish_k
IndStats=array(dim=c(Nrg,3,3))  #3 traits (length,weight,cf) and 3 statistics (mean, cv, n)
sumout=matrix(nrow=Nrg,ncol=Ndetlocs)
for(irg in 1:Nrg){
  irecs=which(rgind==irg)
  for(j in 1:Ndetlocs){
    sumout[irg,j]=sum(CH[irecs,j])
  }

  IndStats[irg,1,1:3]=c(mean(FL[irecs]),sd(FL[irecs])/mean(FL[irecs]),length(irecs))
  IndStats[irg,2,1:3]=c(mean(WGT[irecs]),sd(WGT[irecs])/mean(WGT[irecs]),length(irecs))
  IndStats[irg,3,1:3]=c(mean(CF[irecs]),sd(CF[irecs])/mean(CF[irecs]),length(irecs))
}
write.table(file=here::here("data-raw", "survival_model_STAN", "data.out"),sumout,row.names=RelGp,col.names=c("Release","Woodson","Butte","Sacramento","Delta"))
write.table(file=here::here("data-raw", "survival_model_STAN", "Sz.out"),IndStats,row.names=F,col.names=F)

#Fl=d$fish_length;Wt=d$fish_weight;Cf=d$fish_k #individual covariates that don't vary across reaches
Nseq=25 # # of increments to calculate and plot covariate relationship over

StdRlen=100 #Survival calculated for a standardized reach length of 100 km, then converted to reach specific survival in stan models
Rlen=c(40,88,170,110)
Rmult=Rlen/StdRlen
