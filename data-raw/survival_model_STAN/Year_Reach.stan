//CJS model taken from stan example at https://mc-stan.org/docs/stan-users-guide/mark-recapture-models.html
//This version has same parameter structure as Flora's Year*Reach model
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
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] S_b;//Survival per 10 km in logit space
  array[Nyrs] vector <lower=-10,upper=10> [Nreaches] P_b;
}

transformed parameters{
  array[Nind] vector[Nreaches] surv;
  array[Nind] vector[Ndetlocs] Pcap;
  array[Nind] vector[Ndetlocs] chi; //probability of detecting individual again after first detection by reach
  for(i in 1:Nind){
    for(j in 1:Ndetlocs){  
      if(j<Ndetlocs) surv[i,j]=inv_logit(S_b[yrind[i],j])^Rmult[j];
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
}

generated quantities{
  array[Nyrs] vector[Nreaches] pred_surv;
  array[Nyrs] vector[Nreaches] pred_pcap;
  vector [Nyrs] SurvWoodSac;
  for(iyr in 1:Nyrs){
    for(j in 1:Nreaches){
      pred_surv[iyr,j]=inv_logit(S_b[iyr,j])^Rmult[j];
      pred_pcap[iyr,j]=inv_logit(P_b[iyr,j]); //here pred_pcap[,1] is for woodson, and [,4] is for delta
    }
    SurvWoodSac[iyr]=pred_surv[iyr,2]*pred_surv[iyr,3];
  }
}

