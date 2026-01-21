//CJS model taken from stan example at https://mc-stan.org/docs/stan-users-guide/mark-recapture-models.html
//This version uses year-reach specific pCaps. 
// Survival influenced by fish size at release or 2 or 3-level water year type fixed effect or both cumulatively
//and random effect for year-release group
data {
  int Nind;
  int Nreaches;
  int Ndetlocs;
  int Nyrs;
  int Nrg;
  int NS_bCov; // dimension of discrete covariate
  int UseSizeEffect; // turn on and off the fish size effect
  int Ntribs;//# of tributaries where survival is modelled (Butte + Feather = 2)
  int Nrst; //# of RST considered in mainstem
  int NrstT; //# of RST considered in tribs
  int NindT;// # of individuals released from tributaries summed across Butte and Feather
  int NreachesT;
  int NrgT;// # of release groups summed across Butte and Feather
  int Nsz; // # of size classes
  int NsX; // # of continuous variable such as shasta storage values
  real RmultSac; // distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge
  real RmultTis; // # distance for prediction model from Tisdale to Sacramento per 100km
  vector[NindT] RmultT;//The stream length multiplier for survival (based on distance from tributary release to sac station) for each individual 
  vector[Ntribs] RmultTrib;//by tributary using Boyds release to index the Feather
  array[Nind,Nreaches] real Rmult; // changed this to work for new Rmult
  vector[Nrst] Rmultrst;
  vector[NrstT] RmultTrst;
  array[Nind,Ndetlocs] int CH;
  array[Nind] int yrind;
  array[Nind] int rgind;
  array[Nreaches] int rch_covind;
  array[Nind] int firstCap;
  array[Nind] int lastCap;
  array[Nrg] int rgwy_ind;
  array[NindT] int trib_ind; //The tributary index for each fish
  array[NindT] int firstCapT;// the first detection station of each tributary fish (including release, thus firtCapT = 1 in all cases)
  array[NindT] int lastCapT;// the last detection station for each tributary fish. Note for tributary fish this is 3 or less (Delta = 3 for trib fish)
  array[NindT] int yrindT;// The annual index for each tributary fish for mainstem station pCap mapping
  array[NindT] int rgindT;//The tributary release group index for each tributary fish
  array[NindT,3] int CHT;//The capture history for each tributary fish (digit 1 = release, digit 2 = sac station, digit 3 = delta stations)
  array[NrgT] int trib_rg;//the tributary index for each tributary release group
  array[NrgT] int rgwy_indT;
  array[Nind] real Sz;
  vector[NindT] SzT;
  array[Nsz] real Xsz;
  array[Nind] real CovX; 
  array[NindT] real CovXT;
  array[Nind] real X;
  array[NindT] real XT;
  array[NsX] real Xvec;
  array[NsX] real XvecT;
  
  int Nobs;// # number of individuals*observed detection locations, which allows travel time from release to detection location
  vector[Nobs] ObsTT;
  array[Nind,Nreaches] int TTind; // Index which specifies the element of ObsTT given the 'i' individual in the jth reach.
  int NobsT;// # number of individuals*observed detection locations, which allows travel time from release to detection location
  vector[NobsT] ObsTTT;
  array[NindT,NreachesT] int TTindT; // Index which specifies the element of ObsTT given the 'i' individual in the jth reach.
  vector[Nreaches] ReachKM;
  array[Nind,Nreaches] real ReachKM_ind;
  array[Ntribs,NreachesT] int ReachKMT;
  array[NindT,NreachesT] real ReachKMT_ind;
  vector[Nrst] ReachKMrst;
  vector[NrstT] ReachKMTrst;
  array[NrstT] int TribForRST;
//  array[NrgT,NreachesT] int ReachKMT_rg;
}

parameters {
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] P_b;// reach -specific detection probability
  vector <lower=-10,upper=10> [Nreaches] muPb; // mean detection probability
  vector <lower=0.001> [Nreaches] sdPb; // standard deviation of detection probability
  vector <lower=-5,upper=5> [Nreaches] S_bReach; //survival by reach under baseline covariate
  real <lower=-10,upper=10> S_bCov; //sac mainstem environmental covariate
  real <lower=-10,upper=10> TT_bCov; //sac mainstem environmental covariate for travel time
  array[Nrg, Nreaches] real <lower=-10, upper=10> S_RE; 
  vector <lower=0.001> [Nreaches] RE_sd; //sac mainstem SD of random effect
// vector <lower=-10,upper=10> [Nrg] S_RE; //sac mainstem random year-release group effect
//real <lower=0.001> RE_sd; //sac mainstem SD of random effect
  real <lower=-10,upper=10> S_bSz; // fish size covariate
  real <lower=-10,upper=10> T_bSz; // fish size covariate
  vector <lower=-10,upper=10> [Ntribs] S_bTrib;
  array[NrgT, NreachesT] real  <lower=-10, upper=10> S_REt;
  vector <lower=0.001> [NreachesT] RE_sdT; //tribSD of random effect
 // vector <lower=-5,upper=5> [NrgT] S_REt; //trib random year-release group effect
 // real <lower=0.001> RE_sdT;//tribSD of random effect
  real <lower=-10,upper=10> S_bCovT;//trib environmental covariate
  real <lower=-10,upper=10> TT_bCovT;//trib environmental covariate
   
  real <lower=-10,upper=5> lgB0;//Reach-specific means of travel days per 100 km in log space
  real <lower=0.001, upper=5> Pro_sd;
  array[Nrg, Nreaches] real <lower=-10, upper=10> TT_RE;
   //vector <lower=-10,upper=10> [Nrg] TT_RE;
  vector <lower=0.001> [Nreaches] TTRE_sd;
  //real <lower=0.001> TTRE_sd;
  real <lower=-10,upper=5> lgB0T;//Reach-specific means of travel days per 100 km in log space for tributaries
  real <lower=0.001, upper=5> Pro_sdT;
  array[NrgT, NreachesT] real  <lower=-10, upper=10> TT_RET;
 // vector <lower=-10,upper=10> [NrgT] TT_RET;
 vector <lower=0.001> [NreachesT] TTRE_sdT;
 // real <lower=0.001> TTRE_sdT;
   
  // vector <lower=-10,upper=10> [Nind] S_REind; //sac mainstem random individual effect
  // real <lower=0.001> REind_sd; //sac mainstem SD of random effect
  // vector <lower=-10,upper=10> [NindT] S_REindT; //trib random individual effect
  // real <lower=0.001> REind_sdT; //sac mainstem SD of random effect
}

transformed parameters{
  array[Nind] vector[Nreaches] surv; //survival by reach
  array[Nind] vector[Ndetlocs] Pcap; //detection probability by detection location
  array[Nind] vector[Ndetlocs] chi; //probability of NOT detecting individual again by station
  array[NindT] vector[2] survT; // two reaches release in trib - sac, sac-delta
  array[NindT] vector[3] PcapT;//3 locations: release (not calculated as assumed = 1, 2 = sac, 3 = delta)
  array[NindT] vector[3] chiT; //probability of not detecting individual by station
  
  array[Nind] vector[Nreaches] pTT; //Cummulative travel time from release to detection location
  vector[Nobs] lg_pTT;
  array[NindT] vector[NreachesT] pTTT;  //Cummulative travel time from release to detection location
  vector[NobsT] lg_pTTT;
  real TT;
  
  // Sacramento River fish
  for(i in 1:Nind){
    TT=exp(lgB0 + T_bSz*Sz[i] + TT_bCov*CovX[i] + TT_RE[rgind[i],1]);
    pTT[i,1]=TT*ReachKM_ind[i,1]/100;//predicted travel time is just days/km * km from release to Woodson
    if(TTind[i,1]>0) lg_pTT[TTind[i,1]]=log(pTT[i,1]);//map log predictions to vector for lognormal likelihood calc

    for(j in 2:Nreaches){
      TT=exp(lgB0 + T_bSz*Sz[i] + TT_bCov*CovX[i] + TT_RE[rgind[i],j]);
      pTT[i,j] = pTT[i,j-1] + TT*ReachKM_ind[i,j]/100;//accumulate travel times in downstream direction
      if(TTind[i,j]>0)lg_pTT[TTind[i,j]]=log(pTT[i,j]);
    }
  }
  
  for(i in 1:Nind){
    for(j in 1:Ndetlocs){  
      if(j<Ndetlocs) {
    //   surv[i,j]=inv_logit(S_bReach[rch_covind[j]] + S_bCov*CovX[i] + UseSizeEffect*S_bSz*Sz[i] + S_RE[rgind[i]]+S_REind[i])^Rmult[i,j]; 
      surv[i,j]=inv_logit(S_bReach[rch_covind[j]] + S_bCov*CovX[i] + UseSizeEffect*S_bSz*Sz[i] + S_RE[rgind[i],j])^Rmult[i,j]; 
     }
      if(j>1) Pcap[i,j]=inv_logit(P_b[yrind[i],j-1]);
    }
    chi[i,Ndetlocs]=1.0;
    for(j in 1:Nreaches){
      int r_curr; int r_next;
      r_curr=Ndetlocs-j;  r_next=r_curr+1;
      //            not surviving   or  surviving but not decteced        Accumulate across reaches
      chi[i,r_curr]=(1-surv[i,r_curr]) + surv[i,r_curr]*(1-Pcap[i,r_next]) * chi[i,r_next];        
    }
  }//Nind loop
  
  //Tributary fish (Butte and Feather)
  real TTT;
  for(i in 1:NindT){
    TTT = exp(lgB0T + T_bSz*SzT[i] + TT_bCovT *CovXT[i] + TT_RET[rgindT[i],1]);
    pTTT[i,1]=TTT*ReachKMT_ind[i,1]/100;//predicted travel time is just days/km * km from release to Sacramento
    if(TTindT[i,1]>0) lg_pTTT[TTindT[i,1]]=log(pTTT[i,1]);//map log predictions to vector for lognormal likelihood calc

  for(j in 2:NreachesT){
      TTT=exp(lgB0T + T_bSz*SzT[i] + TT_bCovT*CovXT[i] + TT_RET[rgindT[i],j]);
      pTTT[i,j] = pTTT[i,j-1] + TTT*ReachKMT_ind[i,j]/100;//accumulate travel times in downstream direction
      if(TTindT[i,j]>0)lg_pTTT[TTindT[i,j]]=log(pTTT[i,j]);
    }
}

  for(i in 1:NindT){
    //survival rate in tributary until SAC detection location
     survT[i,1]=inv_logit(S_bTrib[trib_ind[i]] + S_bCovT*CovXT[i] + UseSizeEffect*S_bSz*SzT[i]+ S_REt[rgindT[i],1])^RmultT[i];//RmultT is distance from release point to Sac station
    //   survT[i,1]=inv_logit(S_bTrib[trib_ind[i]] + S_bCovT*CovXT[i] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i]]+S_REindT[i])^RmultT[i];//RmultT is distance from release point to Sac station

    //survival rate in Sac-Delta reach of mainstem
    int j=Nreaches;
    survT[i,2]=inv_logit(S_bReach[rch_covind[j]] + S_bCov*CovX[i] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i],2])^Rmult[1,j];//Same fixed effect for mainstem and tributary fish in mainstem
    //survT[i,2]=inv_logit(S_bReach[rch_covind[j]] + S_bCov*CovX[i] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i]]+S_REindT[i])^Rmult[1,j];//Same fixed effect for mainstem and tributary fish in mainstem
    
    j=4;//pCap at Sac station for tributary fish. Same as for mainstem fish. yrindT assigns the correct year for the trib fish 'i'
    PcapT[i,2]=inv_logit(P_b[yrindT[i],j-1]);//here j points to the mainstem station index, while 2 refers to the digit in the CHT sequence (Sac)
    j=5;//pCap at Delta station for tributary fish
    PcapT[i,3]=inv_logit(P_b[yrindT[i],j-1]);
    
    //Cummulative probability of not being detected by station for tributary fish
    chiT[i,3]=1.0;//probability of not being detected after delta station
    //probability of not being detected after sac station = not surviving from Sac-Delta or suriving Sac-Delta but not detected at Delta station
    chiT[i,2]=(1-survT[i,2]) + survT[i,2]*(1-PcapT[i,3])*chiT[i,3];
    //probability of not being detected after release = not surviving from release in trib to sac station or surviving release-sac but not detected at sac station * probability for next station computed above
    chiT[i,1]=(1-survT[i,1]) + survT[i,1]*(1-PcapT[i,2])*chiT[i,2];
  }
}

model {
    ObsTT[]~lognormal(lg_pTT[],Pro_sd);
   for(j in 1:Nreaches){TT_RE[,j]~normal(0,TTRE_sd[j]);}

   ObsTTT[]~lognormal(lg_pTTT[],Pro_sdT);
    for(j in 1:NreachesT){ TT_RET[,j] ~normal(0,TTRE_sdT[j]);}
    
  for(i in 1:Nind){
    for(j in (firstCap[i]+1):lastCap[i]){  //Loop through all detection stations after first to last station individual was detected at
    1~bernoulli(surv[i,j-1]);//had to be alive prior to last detection station
    CH[i,j]~bernoulli(Pcap[i,j]);
    }
    1~bernoulli(chi[i,lastCap[i]]);//probability of individual never beeing seen again after last detection
  }
  

  for (j in 1:Nreaches){
      S_RE[,j] ~ normal(0,RE_sd[j]);
    }
 // S_RE~normal(0,RE_sd); //RG_sd~gamma(1,1)
//  S_REind~normal(0,REind_sd);
  
  for(j in 1:Nreaches){
    P_b[1:Nyrs,j]~normal(muPb[j],sdPb[j]);
    }
  
  for(i in 1:NindT){
    for(j in (firstCapT[i]+1):lastCapT[i]){//firstCapT is always 1, so +1 means loop starts at 2, and can end at 2 (if lastCapT=2 = Sac) or 3 (delta) 
    1~bernoulli(survT[i,j-1]);//had to be alive prior to last detection station
    CHT[i,j]~bernoulli(PcapT[i,j]);//note j can be 2 or 3
    }
    1~bernoulli(chiT[i,lastCapT[i]]);//probability of individual never beeing seen again after last detection
  }
  
  for (j in 1:NreachesT){
        S_REt[,j] ~ normal(0,RE_sdT[j]);
    }
// S_REt~normal(0,RE_sdT);//random effects for tributary survival for each tributary individual
// S_REindT~normal(0,REind_sdT);
}


generated quantities{

  vector[Nind + NindT] log_lik;

  for(i in 1:Nind){
    log_lik[i]=0;
    for(j in (firstCap[i]+1):lastCap[i]){  //Loop through all detection stations after first to last station individual was detected at
    log_lik[i] += log(surv[i,j-1]);//had to be alive prior to last detection station
    if (CH[i,j]==1){
      log_lik[i] += log(Pcap[i,j]);
    } else {
      log_lik[i] += log1m(Pcap[i,j]);
    }
    }
    log_lik[i] += log(chi[i,lastCap[i]]);//probability of individual never beeing seen again after last detection
  }

  int k = Nind;
  for(i in 1:NindT){
        k=k+1;
        log_lik[k]=0;
    for(j in (firstCapT[i]+1):lastCapT[i]){  //Loop through all detection stations after first to last station individual was detected at
    log_lik[k] += log(survT[i,j-1]);//had to be alive prior to last detection station
    if (CHT[i,j]==1){
      log_lik[k] += log(PcapT[i,j]);
    } else {
      log_lik[k] += log1m(PcapT[i,j]);
    }
    }
     log_lik[k] += log(chiT[i,lastCapT[i]]);//probability of individual never beeing seen again after last detection
  }
  

  array[Nyrs] vector[Nreaches] pred_pcap;
  for(iyr in 1:Nyrs){
    for(j in 2:Ndetlocs){
      pred_pcap[iyr,j-1]=inv_logit(P_b[iyr,j-1]);  
    }
  }
  
  array[Nind] vector[Nreaches] pred_surv;//survival rate by individual
  array[Nind] vector[Nreaches] pred_surv_per100;
 // vector[Nind] SurvWoodSac;
 vector[Nind] SurvRelSac;
 real B0;
 B0=exp(lgB0); //Travel days per 100 km. Needed for output only
 vector[Nind] TT_RelSac;
 array[Nind] vector[Nreaches] TT_reach;
  for(i in 1:Nind){
    TT_reach[i,1]= exp(lgB0 + TT_bCov*X[i] + TT_RE[rgind[i],1])*ReachKM_ind[i,1]/100;
    TT_RelSac[i]= exp(lgB0 + TT_bCov*X[i] + TT_RE[rgind[i],1])*ReachKM_ind[i,1]/100*
                  exp(lgB0 + TT_bCov*X[i] + TT_RE[rgind[i],2])*ReachKM_ind[i,2]/100*
                  exp(lgB0 + TT_bCov*X[i] + TT_RE[rgind[i],3])*ReachKM_ind[i,3]/100;
    pred_surv[i,1]=inv_logit(S_bReach[rch_covind[1]] + S_bCov*X[i]+  S_RE[rgind[i],1])^RmultSac;
    pred_surv_per100[i,1]=inv_logit(S_bReach[rch_covind[1]] +  S_bCov*X[i] + S_RE[rgind[i],1]);

    for(j in 2:Nreaches){ 
     pred_surv[i,j]=inv_logit(S_bReach[rch_covind[j]] +  S_bCov*X[i] + S_RE[rgind[i],j])^Rmult[1,j];
     TT_reach[i,j]= exp(lgB0 + TT_bCov*X[i] + TT_RE[rgind[i],j])*ReachKM_ind[i,j]/100;
     pred_surv_per100[i,j]=inv_logit(S_bReach[rch_covind[j]] +  S_bCov*X[i] + S_RE[rgind[i],j]);
    }
 //   SurvWoodSac[i]=pred_surv[i,2]*pred_surv[i,3]; 
    SurvRelSac[i] = pred_surv[i,1]*pred_surv[i,2]*pred_surv[i,3];
  }
  
//  array[Nsz,NsX] real SurvWoodSacSz;//woodson-sac survival rate by size class and covariate value
  array[Nsz,NsX] real SurvRelSacSz;//release-sac survival rate by size class and covariate value
  

  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
      for(j in 1:NsX){
 //       SurvWoodSacSz[ix,j]=inv_logit(S_bReach[2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,2] * 
//                            inv_logit(S_bReach[3] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,3]; 
        SurvRelSacSz[ix,j]= inv_logit(S_bReach[1] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] )^Rmult[1,1] *  
                              inv_logit(S_bReach[2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,2] * 
                              inv_logit(S_bReach[3] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,3];  
      }
    }
  }
  
  vector[NsX] SurvForecast;
  vector[NsX] SurvForecast_nore;
  vector[NsX] TTForecast;
//  array [Nrst]vector[NsX] TTForecast_rst;
//  array [Nrst] vector [NsX] SurvForecast_rst;
  real s_re1=normal_rng(0,RE_sd[1]); // this can be limited to a max value
  real s_re2=normal_rng(0,RE_sd[2]);
  real s_re3=normal_rng(0,RE_sd[3]);
  real tt_re1=normal_rng(0,TTRE_sd[1]);
  real tt_re2=normal_rng(0,TTRE_sd[2]);
  real tt_re3=normal_rng(0,TTRE_sd[3]);
  for(j in 1:NsX){
    SurvForecast[j]= inv_logit(S_bReach[1] + S_bCov*Xvec[j] + s_re1)^RmultSac* 
                     inv_logit(S_bReach[2] + S_bCov*Xvec[j] + s_re2)^Rmult[1,2] *
                      inv_logit(S_bReach[3] + S_bCov*Xvec[j] + s_re3)^Rmult[1,3];
    SurvForecast_nore[j]= inv_logit(S_bReach[1] + S_bCov*Xvec[j])^RmultSac* 
                     inv_logit(S_bReach[2] + S_bCov*Xvec[j])^Rmult[1,2] *
                      inv_logit(S_bReach[3] + S_bCov*Xvec[j])^Rmult[1,3];
    TTForecast[j] = exp(lgB0 + TT_bCov *Xvec[j] + tt_re1)*ReachKM[1]/100 *
                    exp(lgB0 + TT_bCov *Xvec[j] + tt_re2)*ReachKM[2]/100*
                    exp(lgB0 + TT_bCov *Xvec[j] + tt_re3)*ReachKM[3]/100;
//    for (i in 1:Nrst){
//      SurvForecast_rst[i,j]= inv_logit(S_bReach[1] + S_bCov*Xvec[j] + s_re1)^Rmultrst[i] * 
//                     inv_logit(S_bReach[2] + S_bCov*Xvec[j] + s_re2)^Rmult[1,2] *
//                      inv_logit(S_bReach[3] + S_bCov*Xvec[j] + s_re3)^Rmult[1,3];

//      TTForecast_rst[i,j] = exp(lgB0 + TT_bCov *Xvec[j] + tt_re1)*ReachKMrst[i]/100 *
//                    exp(lgB0 + TT_bCov *Xvec[j] + tt_re2)*ReachKM[2]/100*
//                    exp(lgB0 + TT_bCov *Xvec[j] + tt_re3)*ReachKM[3]/100;
//    }
  }
  
  array[Nsz] vector[NsX] SurvForecastSz; 
  array[Nsz] vector[NsX] SurvForecastSz_nore; 
  array[Nsz] vector[NsX] TTForecastSz;
  array[Nsz,NsX,Nrst] real TTForecastSz_rst;
  array[Nsz,NsX,Nrst] real SurvForecastSz_rst; 
  if(UseSizeEffect == 1){
      for(ix in 1:Nsz){
        for(j in 1:NsX){
          SurvForecastSz[ix,j]= inv_logit(S_bReach[1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re1)^RmultSac * 
                              inv_logit(S_bReach[2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
                              inv_logit(S_bReach[3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re3)^Rmult[1,3];
          SurvForecastSz_nore[ix,j]= inv_logit(S_bReach[1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix])^RmultSac * 
                              inv_logit(S_bReach[2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,2] *
                              inv_logit(S_bReach[3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix])^Rmult[1,3];
          TTForecastSz[ix,j]= exp(lgB0+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re1)*ReachKM[1]/100 *
                           exp(lgB0+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re2)*ReachKM[2]/100 *
                           exp(lgB0+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re3)*ReachKM[3]/100;
            for(i in 1:Nrst){
               SurvForecastSz_rst[ix,j,i]= inv_logit(S_bReach[1] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re1)^Rmultrst[i] *
                              inv_logit(S_bReach[2] + S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
                              inv_logit(S_bReach[3] +  S_bCov*Xvec[j] + S_bSz*Xsz[ix] + s_re3)^Rmult[1,3];
              TTForecastSz_rst[ix,j,i]= exp(lgB0+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re1)*ReachKMrst[i]/100 *
                           exp(lgB0+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re2)*ReachKM[2]/100 *
                           exp(lgB0+TT_bCov*Xvec[j] + T_bSz*Xsz[ix] + tt_re3)*ReachKM[3]/100;
        }
      }
    }
  }
  
array[NindT] vector[2] pred_survT;//Survival of Butte and Feather fish in tributary and sac-delta reach
array[NindT] vector[2] pred_survT_per100;
real B0T;
 B0T=exp(lgB0T); //Travel days per 100 km. Needed for output only
// vector[NindT] TT_RelSacT;
 array[NindT] vector[2] TT_reachT;
  for(i in 1:NindT){
     pred_survT[i,1]=inv_logit(S_bTrib[trib_ind[i]] + S_bCovT*XT[i] + S_REt[rgindT[i],1])^RmultTrib[trib_ind[i]];//survival rate from release to Sac station
     pred_survT[i,2]=inv_logit(S_bReach[rch_covind[Nreaches]] + S_bCov*XT[i] + S_REt[rgindT[i],2] )^Rmult[1,Nreaches];//survival rate from Sac-Delta stations
    
     pred_survT_per100[i,1]=inv_logit(S_bTrib[trib_ind[i]] + S_bCovT*XT[i] + S_REt[rgindT[i],1]);//survival rate from release to Sac station
     pred_survT_per100[i,2]=inv_logit(S_bReach[rch_covind[Nreaches]] + S_bCov*XT[i] + S_REt[rgindT[i],2] );//survival rate from Sac-Delta stations

//     TT_RelSacT[i]= exp(lgB0T + TT_bCovT *XT[i] +TT_RET[rgindT[i]])*sum(ReachKMT_ind[i,1:2])/100;
     TT_reachT[i,1]= exp(lgB0T + TT_bCovT *XT[i] +TT_RET[rgindT[i],1])*ReachKMT_ind[i,1]/100;
     TT_reachT[i,2]= exp(lgB0T + TT_bCovT *XT[i] +TT_RET[rgindT[i],2])*ReachKMT_ind[i,2]/100;
  }
  
  array[Nsz,NsX,Ntribs] real pred_survTSz;//trib release-sac survival rate by size class and covariate value

  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
      for(itrib in 1:Ntribs){
        for(j in 1:NsX){
          pred_survTSz[ix,j,itrib]=inv_logit(S_bTrib[itrib] + S_bCovT*XvecT[j] + S_bSz*Xsz[ix])^RmultTrib[itrib]; // check X term: should it be Xvec? Also missing S_REt?
        }
      }
    }
  }
  
  
array[NsX] vector[Ntribs] TribSurvForecast;//Forecasted survival from Tributary release to Sac station (or Delta if last bit uncommented
array[NsX] vector[Ntribs] TribSurvForecast_nore;
array[NsX] vector[Ntribs] TTForecastT;
//array[NrstT,Ntribs,NsX] TribSurvForecast_rst;//Forecasted survival from Tributary release to Sac station (or Delta if last bit uncommented
//array[NrstT] vector[NsX] TTForecastT_rst;
real s_ret=normal_rng(0,RE_sdT[1]);
real tt_reT=normal_rng(0,TTRE_sdT[1]);
for(itrib in 1:Ntribs){
  for(j in 1:NsX){
        TribSurvForecast[j,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT*XvecT[j] + s_ret)^RmultTrib[itrib]); 
        TribSurvForecast_nore[j,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT*XvecT[j] )^RmultTrib[itrib]); 
        TTForecastT[j,itrib] = exp(lgB0T + TT_bCovT *XvecT[j] + tt_reT)*ReachKMT[itrib,1]/100;
    }
  }
  
  array[Nsz,NsX,Ntribs] real TribSurvForecastSz;
  array[Nsz,NsX,Ntribs] real TribSurvForecastSz_nore;
  array[Nsz,NsX,Ntribs] real TTForecastTSz;
  array[Nsz,NsX,NrstT] real TribSurvForecastSz_rst;
  array[Nsz,NsX,NrstT] real TTForecastTSz_rst;
  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
        for(j in 1:NsX){
            for(itrib in 1:Ntribs){
               TribSurvForecastSz[ix,j,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT*XvecT[j] + S_bSz*Xsz[ix]+ s_ret)^RmultTrib[itrib]);  
               TribSurvForecastSz_nore[ix,j,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT*XvecT[j] + S_bSz*Xsz[ix])^RmultTrib[itrib]);  
               TTForecastTSz[ix,j,itrib] = exp(lgB0T + TT_bCovT *XvecT[j] + T_bSz*Xsz[ix]+ tt_reT)*ReachKMT[itrib,1]/100;
            }
               for(i in 1:NrstT){
                  TTForecastTSz_rst[ix,j,i] = exp(lgB0T + TT_bCovT *XvecT[j] + T_bSz*Xsz[ix]+ tt_reT)*ReachKMTrst[i]/100;
                  TribSurvForecastSz_rst[ix,j,i]=(inv_logit(S_bTrib[TribForRST[i]] + S_bCovT*XvecT[j] + S_bSz*Xsz[ix]+ s_ret)^RmultTrst[i]);
          }
        }
      }
    }
  }

