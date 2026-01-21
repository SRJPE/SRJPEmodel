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
  int NreachesT;
  int NindT;// # of individuals released from tributaries summed across Butte and Feather
  int NrgT;// # of release groups summed across Butte and Feather
  int Nsz; // # of size classes
  int NsX; // # of continuous variable such as shasta storage values
  real RmultSac; // distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge
  real RmultTis; // # distance for prediction model from Tisdale to Sacramento per 100km
  vector[NindT] RmultT;//The stream length multiplier for survival (based on distance from tributary release to sac station) for each individual 
  vector[Ntribs] RmultTrib;//by tributary using Boyds release to index the Feather
  array[Nind,Nreaches] real Rmult; // changed this to work for new Rmult
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
  array[Nind] int CovX; 
  array[NindT] int CovXT;
  array[Nrg] real X;
  array[NrgT] real XT;
  array[NsX] real Xvec;
  array[NsX] real XvecT;

}

parameters {
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] P_b;// reach -specific detection probability
  vector <lower=-10,upper=10> [Nreaches] muPb; // mean detection probability
  vector <lower=0.001> [Nreaches] sdPb; // standard deviation of detection probability
  vector <lower=-10,upper=10> [Nreaches] S_bReach; //survival by reach under baseline covariate
  vector <lower=-10,upper=10> [NS_bCov] S_bCov; //sac mainstem environmental covariate
  array[Nrg, Nreaches] real <lower=-10, upper=10> S_RE; 
  vector <lower=0.001> [Nreaches] RE_sd; //sac mainstem SD of random effect
  real <lower=-10,upper=10> S_bSz; // fish size covariate
  vector <lower=-10,upper=10> [Ntribs] S_bTrib;
  array[NrgT, NreachesT] real  <lower=-10, upper=10> S_REt;
  vector <lower=0.001> [NreachesT] RE_sdT; //tribSD of random effect
  vector<lower=-10,upper=10> [NS_bCov] S_bCovT;//trib environmental covariate
}

transformed parameters{
  array[Nind] vector[Nreaches] surv; //survival by reach
  array[Nind] vector[Ndetlocs] Pcap; //detection probability by detection location
  array[Nind] vector[Ndetlocs] chi; //probability of NOT detecting individual again by station
  array[NindT] vector[2] survT; // two reaches release in trib - sac, sac-delta
  array[NindT] vector[3] PcapT;//3 locations: release (not calculated as assumed = 1, 2 = sac, 3 = delta)
  array[NindT] vector[3] chiT; //probability of not detecting individual by station
  
  // Sacramento River fish
  for(i in 1:Nind){
    for(j in 1:Ndetlocs){  
      if(j<Ndetlocs) {
        if(CovX[i]==0){
          surv[i,j]=inv_logit(S_bReach[rch_covind[j]] + UseSizeEffect*S_bSz*Sz[i] + S_RE[rgind[i],j])^Rmult[i,j]; 
        } else { 
           surv[i,j]=inv_logit(S_bReach[rch_covind[j]] + S_bCov[CovX[i]] + UseSizeEffect*S_bSz*Sz[i] + S_RE[rgind[i],j])^Rmult[i,j];
        }
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
  for(i in 1:NindT){
    //survival rate in tributary until SAC detection location     
    if(CovXT[i]==0){
      survT[i,1]=inv_logit(S_bTrib[trib_ind[i]] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i],1])^RmultT[i];
    } else {
      survT[i,1]=inv_logit(S_bTrib[trib_ind[i]] + S_bCovT[CovXT[i]] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i],1])^RmultT[i];//RmultT is distance from release point to Sac station
    }
    //survival rate in Sac-Delta reach of mainstem
    int j=Nreaches;
    if(CovXT[i]==0){
      survT[i,2]=inv_logit(S_bReach[rch_covind[j]] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i],2])^Rmult[1,j];//
    }else {
      //survT[i,2]=inv_logit(S_bReach[rch_covind[j]] + S_bCov[CovX[i]] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i]])^Rmult[1,j];//Same fixed effect for mainstem and tributary fish in mainstem 
     survT[i,2]=inv_logit(S_bReach[rch_covind[j]] + S_bCov[CovXT[i]] + UseSizeEffect*S_bSz*SzT[i] + S_REt[rgindT[i],2])^Rmult[1,j];//Same fixed effect for mainstem and tributary fish in mainstem 

    }
    
    
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
    
  for(j in 1:Nreaches){P_b[1:Nyrs,j]~normal(muPb[j],sdPb[j]);}
  
  for(i in 1:NindT){
    for(j in (firstCapT[i]+1):lastCapT[i]){//firstCapT is always 1, so +1 means loop starts at 2, and can end at 2 (if lastCapT=2 = Sac) or 3 (delta) 
    1~bernoulli(survT[i,j-1]);//had to be alive prior to last detection station
    CHT[i,j]~bernoulli(PcapT[i,j]);//note j can be 2 or 3
    }
    1~bernoulli(chiT[i,lastCapT[i]]);//probability of individual never beeing seen again after last detection
  }

  for (j in 1:NreachesT){
        S_REt[,j] ~ normal(0,RE_sdT[j]);
    } //random effects for tributary survival for each tributary release group
  
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
  
  array[Nrg] vector[Nreaches] pred_surv;//survival rate by release group
  vector[Nrg] SurvRelSac;
    for(irg in 1:Nrg){
       if(rgwy_ind[irg]==0){
          pred_surv[irg,1]=inv_logit(S_bReach[rch_covind[1]] +  S_RE[irg,1])^RmultSac;
       } else {
         pred_surv[irg,1]=inv_logit(S_bReach[rch_covind[1]] +  S_bCov[rgwy_ind[irg]] +  S_RE[irg,1])^RmultSac;
       }
       
    for(j in 2:Nreaches){ # was 1 before
    if(rgwy_ind[irg]==0){
      pred_surv[irg,j]=inv_logit(S_bReach[rch_covind[j]] +  S_RE[irg,j])^Rmult[1,j];
    } else {
      pred_surv[irg,j]=inv_logit(S_bReach[rch_covind[j]] +  S_bCov[rgwy_ind[irg]] + S_RE[irg,j])^Rmult[1,j];
    }
    }
    SurvRelSac[irg] = pred_surv[irg,1]*pred_surv[irg,2]*pred_surv[irg,3];
  }
  
  array[Nsz,NS_bCov+1] real SurvRelSacSz;//release-sac survival rate by size class and water year type
  
  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
      SurvRelSacSz[ix,1]=inv_logit(S_bReach[1] + S_bSz*Xsz[ix] )^Rmult[1,1] * 
                         inv_logit(S_bReach[2] + S_bSz*Xsz[ix] )^Rmult[1,2] * 
                         inv_logit(S_bReach[3] +  S_bSz*Xsz[ix])^Rmult[1,3];  
      
      for(j in 1:NS_bCov){
        SurvRelSacSz[ix,j+1]= inv_logit(S_bReach[1] + S_bCov[j] + S_bSz*Xsz[ix] )^Rmult[1,1] *  
                              inv_logit(S_bReach[2] + S_bCov[j] + S_bSz*Xsz[ix])^Rmult[1,2] * 
                              inv_logit(S_bReach[3] + S_bCov[j] + S_bSz*Xsz[ix])^Rmult[1,3];  
      }
    }
  }
  
  vector[NS_bCov+1] SurvForecast;
  vector[NS_bCov+1] SurvForecast_nore;
  real s_re1=normal_rng(0,RE_sd[1]);
  real s_re2=normal_rng(0,RE_sd[2]);
  real s_re3=normal_rng(0,RE_sd[3]);
#  vector[NS_bCov+1] SurvTisForecast;
  SurvForecast[1]=  inv_logit(S_bReach[1] + s_re1)^RmultSac * 
                    inv_logit(S_bReach[2] + s_re2)^Rmult[1,2] * 
                    inv_logit(S_bReach[3] + s_re3)^Rmult[1,3];
  
#  SurvTisForecast[1]= inv_logit(S_bReach[3] + s_re)^RmultTis;
 
  for(j in 1:NS_bCov){
    SurvForecast[j+1]= inv_logit(S_bReach[1] + S_bCov[j] + s_re1)^RmultSac * 
                        inv_logit(S_bReach[2] + S_bCov[j] + s_re2)^Rmult[1,2] *
                        inv_logit(S_bReach[3] + S_bCov[j] + s_re3)^Rmult[1,3];
    
    SurvForecast_nore[j+1]= inv_logit(S_bReach[1] + S_bCov[j])^RmultSac * 
                        inv_logit(S_bReach[2] + S_bCov[j])^Rmult[1,2] *
                        inv_logit(S_bReach[3] + S_bCov[j])^Rmult[1,3];
#    SurvTisForecast[j+1]= inv_logit(S_bReach[3] + S_bCov[j] + s_re)^RmultTis;
    
  }
  
  array[Nsz] vector[NS_bCov+1] SurvForecastSz;
  array[Nsz] vector[NS_bCov+1] SurvForecastSz_nore;
#  array[Nsz] vector[NS_bCov+1] SurvTisForecastSz;
  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
      SurvForecastSz[ix,1]= inv_logit(S_bReach[1] + S_bSz*Xsz[ix] + s_re1)^RmultSac * 
                            inv_logit(S_bReach[2] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
                            inv_logit(S_bReach[3] +  S_bSz*Xsz[ix] + s_re3)^Rmult[1,3];
      
      SurvForecastSz_nore[ix,1]= inv_logit(S_bReach[1] + S_bSz*Xsz[ix])^RmultSac * 
                            inv_logit(S_bReach[2] + S_bSz*Xsz[ix])^Rmult[1,2] *
                            inv_logit(S_bReach[3] +  S_bSz*Xsz[ix])^Rmult[1,3];
  #    SurvTisForecastSz[ix,1]= inv_logit(S_bReach[3] +  S_bSz*Xsz[ix] + s_re)^RmultTis;
      
      for(j in 1:NS_bCov){
        SurvForecastSz[ix,j+1]= inv_logit(S_bReach[1] +  S_bCov[j] + S_bSz*Xsz[ix] + s_re1)^RmultSac * 
                                inv_logit(S_bReach[2] + S_bCov[j] + S_bSz*Xsz[ix] + s_re2)^Rmult[1,2] *
                                inv_logit(S_bReach[3] +  S_bCov[j] + S_bSz*Xsz[ix] + s_re3)^Rmult[1,3];
                                
       SurvForecastSz_nore[ix,j+1]= inv_logit(S_bReach[1] +  S_bCov[j] + S_bSz*Xsz[ix])^RmultSac * 
                                inv_logit(S_bReach[2] + S_bCov[j] + S_bSz*Xsz[ix])^Rmult[1,2] *
                                inv_logit(S_bReach[3] +  S_bCov[j] + S_bSz*Xsz[ix])^Rmult[1,3];
        
#        SurvTisForecastSz[ix,j+1]= inv_logit(S_bReach[3] +  S_bSz*Xsz[ix] + s_re)^RmultTis;
        
      }
    }
  }
  
  array[NrgT] vector[2] pred_survT;//Survival of Butte and Feather fish in tributary and sac-delta reach
  for(irg in 1:NrgT){
    if(rgwy_indT[irg]==0){ 
      pred_survT[irg,1]=inv_logit(S_bTrib[trib_rg[irg]] +  S_REt[irg,1])^RmultTrib[trib_rg[irg]];//survival rate from release to Sac station
      pred_survT[irg,2]=inv_logit(S_bReach[rch_covind[Nreaches]] + S_REt[irg,2] )^Rmult[1,Nreaches];//survival rate from Sac-Delta stations
    } else {
      pred_survT[irg,1]=inv_logit(S_bTrib[trib_rg[irg]] + S_bCovT[rgwy_indT[irg]] + S_REt[irg,1])^RmultTrib[trib_rg[irg]];//survival rate from release to Sac station
      pred_survT[irg,2]=inv_logit(S_bReach[rch_covind[Nreaches]] + S_bCov[rgwy_indT[irg]] + S_REt[irg,2] )^Rmult[1,Nreaches];//survival rate from Sac-Delta stations
    }
  }
  
  array[Nsz,NS_bCov+1,Ntribs] real pred_survTSz;//trib release-sac survival rate by size class and water year type

  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
      for(itrib in 1:Ntribs){
        pred_survTSz[ix,1,itrib]=inv_logit(S_bTrib[itrib] + S_bSz*Xsz[ix] )^RmultTrib[itrib];

        for(j in 1:NS_bCov){
          pred_survTSz[ix,j+1,itrib]=inv_logit(S_bTrib[itrib] + S_bCovT[j] + S_bSz*Xsz[ix])^RmultTrib[itrib];
        }
      }
    }
  }
  
  
  array[NS_bCov+1] vector[Ntribs] TribSurvForecast;//Forecasted survival from Tributary release to Sac station (or Delta if last bit uncommented
  array[NS_bCov+1] vector[Ntribs] TribSurvForecast_nore;
  real s_ret=normal_rng(0,RE_sdT[1]);
  for(itrib in 1:Ntribs){
    TribSurvForecast[1,itrib]=(inv_logit(S_bTrib[itrib] + s_ret)^RmultTrib[itrib]); 
        TribSurvForecast[1,itrib]=(inv_logit(S_bTrib[itrib])^RmultTrib[itrib]); 
    for(j in 1:NS_bCov){
      TribSurvForecast[j+1,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT[j] + s_ret)^RmultTrib[itrib]);  
      TribSurvForecast[j+1,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT[j])^RmultTrib[itrib]);  
    }
  }
  
  array[Nsz,NS_bCov+1,Ntribs] real TribSurvForecastSz;
  array[Nsz,NS_bCov+1,Ntribs] real TribSurvForecastSz_nore;
  if(UseSizeEffect == 1){
    for(ix in 1:Nsz){
      for(itrib in 1:Ntribs){
        TribSurvForecastSz[ix,1,itrib]=(inv_logit(S_bTrib[itrib]+ S_bSz*Xsz[ix]+ s_ret)^RmultTrib[itrib]); 
        TribSurvForecastSz_nore[ix,1,itrib]=(inv_logit(S_bTrib[itrib]+ S_bSz*Xsz[ix])^RmultTrib[itrib]); 
        for(j in 1:NS_bCov){
          TribSurvForecastSz[ix,j+1,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT[j] + S_bSz*Xsz[ix]+ s_ret)^RmultTrib[itrib]); 
          TribSurvForecastSz_nore[ix,j+1,itrib]=(inv_logit(S_bTrib[itrib] + S_bCovT[j] + S_bSz*Xsz[ix])^RmultTrib[itrib]);
        }
      }
    }
  }  
}


