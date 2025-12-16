library(rstan)
library(tidyverse)
library(SRJPEmodel)
library(SRJPEdata)
library("R2WinBUGS")

# liz add
setwd(here::here("data-raw", "juvenile_abundance", "testing"))

RunpCapModel=T
RunMainstem=F

d0a=SRJPEdata::years_to_include_rst_data |>
  mutate(exclude = case_when(site == "lbc" ~ TRUE,
                             site == "mill creek" & run_year == 2025 ~ TRUE,
                             TRUE ~ FALSE)) |>
  filter(!exclude) |>
  select(-exclude)

if(RunMainstem==F){
  d0=subset(d0a,site!="knights landing"& site!="tisdale" & site!="red bluff diversion dam")
  OutDir="Output_Trib/"
  pCapfn= paste0(OutDir,"pCap_model.rdata")
  min_pCap=0.005

  if(RunpCapModel==T){#if not doing mainstem can run the pcap model out here
    pCap_inputs <- prepare_pCap_inputs(mainstem = RunMainstem)
    #pCap_inputs <- prepare_pCap_inputs(mainstem = RunMainstem, drop_trib_sites=T, sites_to_drop=c("steep riffle","okie dam"))

    Snm=paste0(OutDir,"Inputs_pCap.rdata")
    save(pCap_inputs,file=Snm)

    pCap <- fit_pCap_model(pCap_inputs)
    save(pCap, file = pCapfn)

  } else {
    load(file=pCapfn)
  }

} else {#is a mainstem run
  #will have to run pCap model insite isite loop
  d0=subset(d0a,site=="knights landing" | site=="tisdale")
  OutDir="Output_Main/"
  min_pCap=0.0005

}

uniqsite=unique(d0$site)
Nsites=length(uniqsite)

for(isite in 1:Nsites){
  d1=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]
  uniqyr=unique(sort(d1$run_year))
  Nyrs=length(uniqyr)

  # liz add
  print(DoSite)

  if(RunMainstem==T){
    pCapfn= paste0(OutDir,DoSite,"_pCap_model.rdata")


    if(RunpCapModel==T){
      pCap_inputs <- prepare_pCap_inputs(mainstem = RunMainstem, mainstem_site=DoSite)
      Snm=paste0(OutDir,DoSite,"_Inputs_pCap.rdata")
      save(pCap_inputs,file=Snm)

      pCap <- fit_pCap_model(pCap_inputs)
      save(pCap, file = pCapfn)

    } else { #mainstem model but don't want to run pCap model and instead just load correct file
      load(file=pCapfn)
    }
  }


  for(iyr in 1:Nyrs){#1:Nyrs
  #for(iyr in 7:7){
    DoYr=uniqyr[iyr]

    # liz add
    print(DoYr)

    ModNm=paste0(uniqsite[isite],"_",uniqyr[iyr])


    abundance_inputs <- prepare_abundance_inputs(site = DoSite,run_year = DoYr, effort_adjust = T, pcap_model_object = pCap)
    Snm=paste0(OutDir,"Inputs_",ModNm,".rdata")
    save(abundance_inputs,file=Snm)

    # Setup data forabundance model
    Nstrata=abundance_inputs$inputs$data$Nstrata
    Nstrata_wc=abundance_inputs$inputs$data$Nstrata_wc
    Uwc_ind=abundance_inputs$inputs$data$Uwc_ind
    u=abundance_inputs$inputs$data$u
    K=abundance_inputs$inputs$data$K
    ZP=abundance_inputs$inputs$data$ZP
    lgN_max=abundance_inputs$inputs$data$lgN_max

    lgN_max=rep(log(0.001*(mean(u,na.rm=T)+1)/min_pCap),Nstrata)
    #lgN_max=rep(20,Nstrata)
    for(j in 1:Nstrata_wc){
        if(is.na(u[j])==F) lgN_max[Uwc_ind[j]]=log(0.001*(u[j]+1)/min_pCap)
    }

    lt_pCap_mu=abundance_inputs$lt_pCap_Us$lt_pCap_mu
    lt_pCap_sd=abundance_inputs$lt_pCap_Us$lt_pCap_sd
    lt_pCap_tau<-1/lt_pCap_sd^2

    print(c(uniqyr[iyr],mean(u)))
    data<-list("Nstrata","Nstrata_wc","u","K","ZP","Uwc_ind","lgN_max","lt_pCap_mu","lt_pCap_tau")
    #data<-list(Nstrata=Nstrata,Nstrata_wc=Nstrata_wc,u=u,K=K,ZP=ZP,Uwc_ind=Uwc_ind,lgN_max=lgN_max,lt_pCap_mu=lt_pCap_mu,lt_pCap_tau=lt_pCap_tau)

    ini_b_sp=rep(1,K)
    ini_lgN=rep(log(0.001*(min(u)+1)/0.025),Nstrata)
    for(i in 1:Nstrata_wc) ini_lgN[Uwc_ind[i]]=log(0.001*(u[i]+1)/0.025)
    inits1<-list(b_sp=ini_b_sp,tau.N=1,tau.Ne=1,lg_N=ini_lgN,lt_pCap_U=lt_pCap_mu)
    inits2<-inits1;inits3<-inits1;inits<-list(inits1,inits2,inits3)

    parameters<-c("lt_pCap_U","pCap_U","N","Ntot","sd.N","sd.Ne","lg_CumN")
    Nmcmc=2000;Nburnin=500;Nthin=2;Nchains=3
    abundance<-bugs(data, inits, parameters, "C:/Users/Liz/Documents/SRJPEmodel/model_files/abundance_model.bug", debug=F,n.chains=Nchains, n.burnin=Nburnin,n.thin=Nthin,n.iter=Nmcmc,codaPkg=F,DIC=T,clearWD=T,bugs.directory="C:/Users/Liz/Documents/SRJPEmodel/data-raw/WinBUGS14")	#run bugs model

    #dp=abundance$sims.list$Ntot#as.data.frame(abundance,pars=Vnm)
    #Ntot_cv=round(100*sd(dp)/mean(dp,digits=0))

    Snm=paste0(OutDir,"N_",ModNm,".rdata")
    save(abundance,file=Snm)

  }
}
