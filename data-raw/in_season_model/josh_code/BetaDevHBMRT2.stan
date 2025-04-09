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
  array[Nyrs,2] real CovX;//covariate
  array[2] real UseCovX;//0 or 1 value used to turn off covariate effect
  
  int Nfor;// number of forecasts to make
  array[Nfor] int For_ewk; //proportional week for each forecast
  real For_Nx_mu;//mean abundance through week for each forecast
  real For_Nx_sd;//sd of abundance
  array[2] real For_CovX;
  int NxError;//0 or 1 flag to not simulate or simulate error in Nx
  
}


parameters {
 array[Nyrs] vector[2] RTpars;//column 1 are annual logit_phi values, col 2 are annual log_lambda va;ies
 array[Nwks,Nyrs] real dev;//deviations (random effects) in proportion by week from beta
 real<lower=0.1> sd_pro;//sd for random effects
 array[Nwks,Nyrs] real log_Nx;//estimated state for abundance through week X
 array[2] real <lower=-2,upper=2> bCov; //fixed effects on on annual phi's and lambda's
 
 //Hyper parameters for MVN
 row_vector[2] muRT; // mean to fit for MVN hyper distribution for logit_phi and log_lambda
 array[2] real <lower=0> sigmaRT; // scales to fit (sd for logit_phi across years, sd for log_lambda across years)
 real<lower=-1,upper=1> rho; // correlation to fit (correlation between annual logit_phi and log_lambda values)
}


transformed parameters{
  array[Nyrs] real alpha;
  array[Nyrs] real beta;
  array[Nyrs] real phi;
  array[Nyrs] real lambda;
  
  array[Nwks,Nyrs] real bp;//predicted proportion of run passing on each week (not cummulative) from beta
  array[Nwks,Nyrs] real pdev;
  array[Nwks,Nyrs] real cp;//cummulative probability of total run passing by end of each week

  cov_matrix[2] vcvRT;//setup variance-covariance matrix for annual phi-lambda MVN
  vcvRT[1,1] = square(sigmaRT[1]);//convert from sd to variance. The diagonals (var) of var-covar matrix MVN needs
  vcvRT[2,2] = square(sigmaRT[2]);
  vcvRT[1,2] = rho * sigmaRT[1] * sigmaRT[2];//covariance off-diagonal
  vcvRT[2,1] = vcvRT[1,2];
 
  real bp0;
  for(iyr in 1:Nyrs){
    
    //map MVN parameters to more obvious ones
    phi[iyr]=inv_logit(RTpars[iyr,1] + UseCovX[1]*bCov[1]*CovX[iyr,1]);//add in fixed effect and transform
    lambda[iyr]=exp(RTpars[iyr,2]+ UseCovX[2]*bCov[2]*CovX[iyr,2]);
    
    //convert mean and count of beta into alpha and beta terms
    alpha[iyr]=lambda[iyr] * phi[iyr];
    beta[iyr]= lambda[iyr] * (1 - phi[iyr]); 
    
    real sumpp=0;
    for(iwk in 1:Nwks){
	      // beta probability given weekly  theta and alpha and beta parameters (won't sum to one yet)
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
 
  array[Nwks,Nyrs] real pNtot;//predicted total abundance

  //vectorized MVN which returns annual values of logit_phi (RTpars[1:Nyrs,1]) and 
  //log_lambda (RTpars[1:Nyrs,2])) given means muRT[1:2] and variance covariance matrix vcvRT
  RTpars~multi_normal(muRT,vcvRT);
  
  for(iyr in 1:Nyrs){
    //sd_pro[iyr]~gamma(pr_sd_pro[1],pr_sd_pro[2]);//if no gamma hyper just use a gamma prior
    
    for(iwk in 1:Nestwks[iyr]){
      //estimate true state of Nx in log space (log_Nx) given estimate Nx_mu (in log space) and observation error Nx_sd
      Nx_mu[iwk,iyr]~normal(log_Nx[iwk,iyr],Nx_sd[iwk,iyr]);
      
      if(NxError==0){
        pNtot[iwk,iyr]=exp(Nx_mu[iwk,iyr])/cp[ewk[iwk,iyr],iyr];//no error in estimate of abundance through Wk iwk
      } else {
        pNtot[iwk,iyr]=exp(log_Nx[iwk,iyr])/cp[ewk[iwk,iyr],iyr];//with error in estimate of abundance through wk iwk
      }
      
      //The main data likelihood predicting observed total run size Ntot_mu (in log space) based on mean predicted by model 
      //(pNtot) and observation error in log total run size Ntot_sd 
      Ntot_mu[iyr]~normal(log(pNtot[iwk,iyr]),Ntot_sd[iyr]);
      
    }//iwk
    
    for(iwk in 1:Nwks){dev[iwk,iyr]~normal(0,sd_pro);}//extra beta error in proportion
  }//iyr
  
}

generated quantities{
  real mu_phi=inv_logit(muRT[1]);
  real mu_lambda=exp(muRT[2]);
  
    //define parameters for beta distribution used for forecast
  vector[2] sRT;
  sRT=multi_normal_rng(muRT,vcvRT); //simulate a logit_phi and log_lambda value for the forecast
  real sphi=inv_logit(sRT[1]+UseCovX[1]*bCov[1]*For_CovX[1]);//add in fixed effect and transform
  real slambda=exp(sRT[2]+UseCovX[2]*bCov[2]*For_CovX[2]);
  real salpha=slambda * sphi;//convert to alpha and beta values
  real sbeta=slambda*(1-sphi);
  
  //predict the proportion of run passing by week based on simulated alpha and beta values and simulated deviations
  real sbp0;real sdev;real ssumpp=0;
  array[Nwks] real sbp;
  for(iwk in 1:Nwks){
    sbp0=pow(theta[iwk],salpha-1.0)*pow(1.0-theta[iwk],sbeta-1.0);
    sdev=normal_rng(0,sd_pro);
    sbp[iwk]=sbp0*exp(sdev);
    ssumpp=ssumpp+sbp[iwk];
  }
  
  array[Nwks] real spdev;array[Nwks] real For_cp;
  for(iwk in 1:Nwks){
      spdev[iwk]=sbp[iwk]/ssumpp; //pdf values across weeks will now sum to one
      For_cp[iwk]=sum(spdev[1:iwk]); //create cdf
  }
  
  //Finally, predict total abundance based on abundance through week x and simulated cdf
  array[Nfor] real For_Ntot;
  real slog_Nx;
  if(NxError==0){
    slog_Nx=For_Nx_mu;
  } else {
    slog_Nx=normal_rng(For_Nx_mu,For_Nx_sd);//simulate error in abundance through week X
  }
  
  for(ifor in 1:Nfor){
      For_Ntot[ifor]=exp(slog_Nx)/For_cp[For_ewk[ifor]];
  }
}
