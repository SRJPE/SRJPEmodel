
data {
  int Ntribs;
  int Nfor;
  array[Ntribs,Nfor] real obs_Nx_mu; //mean of log cummulative abundance through each forecast week
  array[Ntribs,Nfor] real obs_Nx_sd; //sd of log cummulative abundance
  vector[Ntribs]  srNtot_mu; //log of mean ofabundance estimate from SR model that will be used as a prior
  vector[Ntribs] srNtot_sd; //log of sd of abundance estimate
  array[Ntribs,Nfor] real cp_mu; // mean of cumulative proportion of run passing by each forecast week in logit space from inseason model
  array[Ntribs,Nfor] real cp_sd; //sd of cumulative proportion
  vector[Ntribs] DS_surv_mu;//mean of weighted RST-Delta survival rate in logit sapace
  vector[Ntribs] DS_surv_sd;//sd
}

parameters {
  array[Ntribs,Nfor] real pred_lgNtot;  //predicted annual abundance in log space
  array[Ntribs,Nfor] real lt_cp;       //predicted cum prop past trap on forecast week in logit space
  array[Ntribs] real<lower=-6.5> lt_DS_surv; //survival to delta in logit space (can't be lower than ~ 0.15% survival)
}

transformed parameters{
   array[Ntribs,Nfor] real pred_lgNx;
   array[Ntribs,Nfor] real cp;
   
   for(itrib in 1:Ntribs){
      for(ifor in 1:Nfor){
        cp[itrib,ifor]=inv_logit(lt_cp[itrib,ifor]);//convert logit cp past trap
        pred_lgNx[itrib,ifor]=log(exp(pred_lgNtot[itrib,ifor])*cp[itrib,ifor]);//predicted abundance on forecast wk is estimated annual abundance * cummulative proportion that has passed through that wk
      }
   }
}

model {//priors and data likelihood
  for(itrib in 1:Ntribs){
    for(ifor in 1:Nfor){
      obs_Nx_mu[itrib,ifor]~normal(pred_lgNx[itrib,ifor],obs_Nx_sd[itrib,ifor]);//'data' likelihood comparing 'observed' and predicted cum abundance through forecast week
      pred_lgNtot[itrib,ifor]~normal(srNtot_mu[itrib],srNtot_sd[itrib]);//prior on annual abundance from SR model
      lt_cp[itrib,ifor]~normal(cp_mu[itrib,ifor],cp_sd[itrib,ifor]);/// prior on cum proportion passing RST by forecast week in logit space from inseason model
    }
    lt_DS_surv[itrib]~normal(DS_surv_mu[itrib],DS_surv_sd[itrib]);
  }
}

generated quantities{
  array[Ntribs,Nfor] real pred_Ntot;
  array[Nfor] real pred_Ntot_all;
  array[Ntribs] real DS_surv;
  array[Ntribs,Nfor] real JPE_trib;
  array[Nfor] real JPE;
  
  for(itrib in 1:Ntribs){
    DS_surv[itrib]=inv_logit(lt_DS_surv[itrib]);//RST to Delta-entry survival rate
  }
  
  for(ifor in 1:Nfor){
    for(itrib in 1:Ntribs){
      pred_Ntot[itrib,ifor]=exp(pred_lgNtot[itrib,ifor]);//convert annual abundance estimate from log space
      JPE_trib[itrib,ifor]=pred_Ntot[itrib,ifor]*DS_surv[itrib];//predict abundance at Delta entry
    }
    pred_Ntot_all[ifor]=sum(pred_Ntot[1:Ntribs,ifor]);//outmigrant abundance at RST summed across tributaries
    JPE[ifor]=sum(JPE_trib[1:Ntribs,ifor]);//outmigrant abundance at delta entry summed across tribs.
  }
}

