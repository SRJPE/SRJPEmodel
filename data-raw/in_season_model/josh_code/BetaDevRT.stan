data {
  int Nyrs;
  int Nwks;
  array[Nyrs] int Nestwks;
  array[Nwks] real theta;
  array[Nwks,Nyrs] int ewk;
  array[Nyrs] real Ntot_mu;//mean of total abundance for run in log space
  array[Nyrs] real Ntot_sd;//sd abundance for run in log space
  array[Nwks,Nyrs] real Nx_mu;//mean of abundance through each week in log space
  array[Nwks,Nyrs] real Nx_sd;//sd of abundance through each weekin log space
  
}


parameters {
 array[Nyrs] real <lower=0, upper=1> phi;   //mean of beta
 //array[Nyrs] real  logit_phi;   //mean of beta
 array[Nyrs] real <lower=0.1> lambda;      // count of beta
 array[Nwks,Nyrs] real dev;
 array[Nyrs] real<lower=0.1> sd_pro;
 array[Nwks,Nyrs] real log_Nx;//estimated state for abundance through month X
 //real mu_phi;
 //real<lower=0.1> sd_logit_phi;
}


transformed parameters{
  array[Nyrs] real alpha;
  array[Nyrs] real beta;
  //array[Nyrs] real phi;
  real mu_logit_phi;
  array[Nwks,Nyrs] real bp;//predicted proportion of run passing on each week (not cummulative) from beta
  array[Nwks,Nyrs] real pdev;
  array[Nwks,Nyrs] real cp;//cummulative probability of total run passing by end of each week

  //mu_logit_phi=logit(mu_phi);
  
  real bp0;
  for(iyr in 1:Nyrs){
    //phi[iyr]=inv_logit(logit_phi[iyr]);
    
    alpha[iyr]=lambda[iyr] * phi[iyr];
    beta[iyr]= lambda[iyr] * (1 - phi[iyr]); 
    
    real sumpp=0;
    for(iwk in 1:Nwks){
	      // beta probability given theta and alpha and beta parameters (won't sume to one yet)
	      bp0=pow(theta[iwk],alpha[iyr]-1.0)*pow(1.0-theta[iwk],beta[iyr]-1.0);
        bp[iwk,iyr]=bp0*exp(dev[iwk,iyr]);//beta with weekly deviations that range from 0 - >>1
        sumpp=sumpp+bp[iwk,iyr];//so bp's sum to 1 based on next loop
      }
      
      for(iwk in 1:Nwks){
        pdev[iwk,iyr]=bp[iwk,iyr]/sumpp; //pdf values across weeks will now sum to one
        cp[iwk,iyr]=sum(pdev[1:iwk,iyr]); //create cdf
      }
  }//iyr
}

model {
 
  array[Nwks,Nyrs] real pNtot;

  for(iyr in 1:Nyrs){
    //logit_phi[iyr]~normal(mu_logit_phi,sd_logit_phi);
    
    for(iwk in 1:Nestwks[iyr]){
      //estimate true state of Nx in log space (log_Nx) given estimate Nx_mu (in log space) and observation error Nx_sd
      Nx_mu[iwk,iyr]~normal(log_Nx[iwk,iyr],Nx_sd[iwk,iyr]);
          
      //The main data likelihood predicting observed total run size Ntot_mu (in log space) based on mean predicted by model 
      //(pNtot) and observation error in log total run size Ntot_sd
      //pNtot[iwk,iyr]=exp(Nx_mu[iwk,iyr])/cp[ewk[iwk,iyr],iyr];//no error in pre-season abundance
      pNtot[iwk,iyr]=exp(log_Nx[iwk,iyr])/cp[ewk[iwk,iyr],iyr];//with error in pre-season abundance
      Ntot_mu[iyr]~normal(log(pNtot[iwk,iyr]),Ntot_sd[iyr]);
      
    }//iwk
    
    //phi[iyr] ~ beta(1, 1); lambda[iyr] ~ pareto(0.1, 1.5); sd_pro[iyr]~gamma(2,5);
    
    for(iwk in 1:Nwks){
      dev[iwk,iyr]~normal(0,sd_pro[iyr]);
    }
  }//iyr
  
}
