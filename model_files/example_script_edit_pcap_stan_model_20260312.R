
# 3-12-2026
# example script for making modifications to STAN pCap model in git workflow
library(tidyverse)
library(SRJPEdata)
#library(SRJPEmodel) # and any others
library(rstan)

options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)
nchains=3;niter=10000

#source("R/juvenile_abundance.R") # get functions

for(ii in 1:3){#1:4
  if(ii==1){
    pCap_inputs <- prepare_pCap_inputs(mainstem=F)
    Mnm="pCap_trib"
    parlist=c("trib_mu_P", "trib_sd_P", "flow_mu_P", "flow_sd_P","logit_pCap", "b0_pCap", "b_flow", "pro_sd_P","yr_sd_P","yr_re","log_lik")
    fnsave=paste0("model_files/Output/",Mnm,".Rdata")

  } else if(ii==2) {
    Mnm="pCap_mainstem_skew_re"
    parlist=c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P","yr_sd_P","yr_re","alpha","log_lik")
    pCap_inputs <- prepare_pCap_inputs(mainstem=T,mainstem_site="knights landing")
    fnsave=paste0("model_files/Output/",Mnm,"_knights landing.Rdata")

  } else if(ii==3) {
    Mnm="pCap_mainstem_skew_re"
    parlist=c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P","yr_sd_P","yr_re","alpha","log_lik")
    pCap_inputs <- prepare_pCap_inputs(mainstem=T,mainstem_site="tisdale")
    fnsave=paste0("model_files/Output/",Mnm,"_tisdale.Rdata")
  }


  pcap <- fit_pCap_model(pCap_inputs) #automates some of arguments, but code below better for testing

  save(pcap,file=fnsave)

}

