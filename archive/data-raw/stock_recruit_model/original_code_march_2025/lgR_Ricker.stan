
data {
  int Nyrs;
  array[Nyrs] real SP;
  array[Nyrs] real mu_obslgR;
  array[Nyrs] real sd_obslgR;
  array[Nyrs] real X;
}

parameters {
 real alpha;
 real <upper=0> beta;
 real <lower=-5,upper=5> gamma;
 real <lower=0.01> sd_pro;
}

transformed parameters{
  array[Nyrs] real pred_R;
  array[Nyrs] real sdTot;
  for(iyr in 1:Nyrs){
    pred_R[iyr] = SP[iyr]*exp(alpha + beta*SP[iyr] + gamma*X[iyr]);
    sdTot[iyr]=pow(pow(sd_obslgR[iyr],2) + pow(sd_pro,2),0.5);
  }
}

model {
  for(iyr in 1:Nyrs){
    mu_obslgR[iyr]~normal(log(pred_R[iyr]),sdTot[iyr]);
  }
  
  alpha~normal(0,1000);
  beta~normal(0,1000);
  gamma~normal(0,1000);
}

