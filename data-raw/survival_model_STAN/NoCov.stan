//CJS model taken from stan example at https://mc-stan.org/docs/stan-users-guide/mark-recapture-models.html
//This version uses year-reach specific pCaps and no covariate effects and random effect for year-release group
data {
  int Nind;
  int Nreaches;
  int Ndetlocs;
  vector[Nreaches] Rmult;
  int Nyrs;
  int Nrg;
  array[Nind,Ndetlocs] int CH;
  int yrind[Nind];
  int rgind [Nind];
  int rch_covind [Nreaches];
  vector[Nind] CovX;
  int firstCap [Nind];
  int lastCap [Nind];
  int rgwy_ind[Nrg];
}

parameters {
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] P_b;
  vector <lower=-10,upper=10> [Nreaches] muPb;
  vector <lower=0.001> [Nreaches] sdPb;
  vector <lower=-10,upper=10> [Nreaches-1] S_bReach;
  vector <lower=-10,upper=10> [Nrg] S_RE;
  real <lower=0.001> RE_sd;
}

transformed parameters{
  array[Nind] vector[Nreaches] surv;
  array[Nind] vector[Ndetlocs] Pcap;
  array[Nind] vector[Ndetlocs] chi; //probability of detecting individual again after first detection by reach
  
  for(i in 1:Nind){
    for(j in 1:Ndetlocs){  
      if(j<Ndetlocs) surv[i,j]=inv_logit(S_bReach[rch_covind[j]] + S_RE[rgind[i]])^Rmult[j];
      if(j>1) Pcap[i,j]=inv_logit(P_b[yrind[i],j-1]);
    }
    
    chi[i,Ndetlocs]=1.0;
    for(j in 1:Nreaches){
      int r_curr; int r_next;
      r_curr=Ndetlocs-j;  r_next=r_curr+1;
      chi[i,r_curr]=(1-surv[i,r_curr]) + surv[i,r_curr]*(1-Pcap[i,r_next]) * chi[i,r_next];         //            not surviving   or  surviving but not decteced     accumulate
    }
  }//Nind loop
}

model {
  for(i in 1:Nind){
      for(j in (firstCap[i]+1):lastCap[i]){  //Loop through all detection stations after first to last station individual was detected at
        1~bernoulli(surv[i,j-1]);//had to be alive prior to last detection station
        CH[i,j]~bernoulli(Pcap[i,j]);
      }
      1~bernoulli(chi[i,lastCap[i]]);//probabily of individual never beeing seen again after being last observed 
  }
  S_RE~normal(0,RE_sd); //RG_sd~gamma(1,1)
  for(j in 1:Nreaches){P_b[1:Nyrs,j]~normal(muPb[j],sdPb[j]);}
}

generated quantities{
  array[Nyrs] vector[Nreaches] pred_pcap;
  for(iyr in 1:Nyrs){
    for(j in 2:Ndetlocs){
      pred_pcap[iyr,j-1]=inv_logit(P_b[iyr,j-1]);  
    }
  }
  
  array[Nrg] vector[Nreaches] pred_surv;
  vector[Nrg] SurvWoodSac;
  for(irg in 1:Nrg){
    for(j in 1:Nreaches){
        pred_surv[irg,j]=inv_logit(S_bReach[rch_covind[j]] +  S_RE[irg])^Rmult[j];
    }
    SurvWoodSac[irg]=pred_surv[irg,2]*pred_surv[irg,3]; 
  }
 
  real SurvForecast;
  real s_re=normal_rng(0,RE_sd);
  SurvForecast=inv_logit(S_bReach[2] +  s_re)^Rmult[2] * inv_logit(S_bReach[3] +  s_re)^Rmult[3]; 
  
}

