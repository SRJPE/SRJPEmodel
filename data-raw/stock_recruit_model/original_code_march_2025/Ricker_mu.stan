
data {
  int Nyrs;
  array[Nyrs] real SP;
  array[Nyrs] real mu_obslgRS;
  array[Nyrs] real sd_obslgRS;
}

parameters {
 real alpha;
 real <lower=0.01> sd_pro;
 array[Nyrs] real dev;
}

transformed parameters{
  array[Nyrs] real lgRS;
  array[Nyrs] real pred_R;
  //array[Nyrs] real sdTot;
  for(iyr in 1:Nyrs){
    lgRS[iyr] = alpha + dev[iyr];
    //sdTot[iyr]=pow(pow(sd_obslgRS[iyr],2) + pow(sd_pro,2),0.5);
    pred_R[iyr]=SP[iyr]*exp(lgRS[iyr]);
  }
}

model {
  for(iyr in 1:Nyrs){
    //mu_obslgRS[iyr]~normal(lgRS[iyr],sdTot[iyr]);
    dev[iyr]~normal(0,sd_pro);
    mu_obslgRS[iyr]~normal(lgRS[iyr],sd_obslgRS[iyr]);
  }
  
  alpha~normal(0,1000);
}

generated quantities{
  vector[Nyrs] log_lik;
  for (iyr in 1:Nyrs) {
    //log_lik[iyr] = normal_lpdf(mu_obslgRS[iyr] | lgRS[iyr],sdTot[iyr]);
    log_lik[iyr] = normal_lpdf(mu_obslgRS[iyr] | lgRS[iyr],sd_obslgRS[iyr]);
  }
  
}

