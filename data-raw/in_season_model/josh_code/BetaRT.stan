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
 array[Nyrs] real <lower=0.1> lambda;      // count of beta
 array[Nwks,Nyrs] real pdev;
 array[Nyrs] real<lower=0.1> sd_pro;
 array[Nwks,Nyrs] real log_Nx;//estimated state for abundance through month X
}


transformed parameters{
  
  array[Nyrs] real alpha;
  array[Nyrs] real beta;
  array[Nwks] real pp;//predicted proportion of run passing on each week (not cummulative)
  array[Nwks,Nyrs] real cp;//cummulative probability of total run passing by end of each week
  
  
  for(iyr in 1:Nyrs){
    alpha[iyr]=lambda[iyr] * phi[iyr];
    beta[iyr]= lambda[iyr] * (1 - phi[iyr]); 
    
    real sumpp=0;
    for(iwk in 1:Nwks){
      pp[iwk]=pow(theta[iwk],alpha[iyr]-1.0)*pow(1.0-theta[iwk],beta[iyr]-1.0);
      sumpp=sumpp+pp[iwk];
      }
    for(iwk in 1:Nwks){
      pp[iwk]=pp[iwk]/sumpp; //pdf values across weeks will now sum to one
      cp[iwk,iyr]=sum(pp[1:iwk]); //create cdd
    }
    
    //Same answer as above but takes241 seconds for 1000 iterations vs. 40 seconds based on code above
	  //for(iwk in 1:Nwks){cp[iwk,iyr]=beta_cdf(theta[iwk] | alpha[iyr], beta[iyr]);}
	  
  }//iyr
}

model {
 
  array[Nwks,Nyrs] real pNtot;

  for(iyr in 1:Nyrs){
    for(iwk in 1:Nestwks[iyr]){
      //estimate true state of Nx in log space (log_Nx) given observation error Nx_sd
      Nx_mu[iwk,iyr]~normal(log_Nx[iwk,iyr],Nx_sd[iwk,iyr]);
          
      //The main data likelyood predicting observed total run size Ntot_mu (in log space) based on mean predicted by model 
      //(pNtot) and observation error in log total run size Ntot_sd
      //pNtot[iwk,iyr]=exp(Nx_mu[iwk,iyr])/cp[ewk[iwk,iyr],iyr];//no error in pre-season abundance
      pNtot[iwk,iyr]=exp(log_Nx[iwk,iyr])/cp[ewk[iwk,iyr],iyr];//with error in pre-season abundance
      Ntot_mu[iyr]~normal(log(pNtot[iwk,iyr]),Ntot_sd[iyr]);
      
    }//iwk
    
    phi[iyr] ~ beta(1, 1); 
    lambda[iyr] ~ pareto(0.1, 1.5);
    
  }//iyr
  
}
