# 3-27-2026
# example script for making modifications to STAN pCap model in git workflow
library(SRJPEdata)
library(SRJPEmodel) # and any others
library(tidyverse)
library(dplyr)
library(brms)
library("R2WinBUGS")


OutDir="C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/output/0p5minP"

d0a=SRJPEdata::years_to_include_rst_data
d0=subset(d0a,site!="lbc" & site!="red bluff diversion dam")

siteO=c("ucc","lcc","ubc","mill creek","deer creek","okie dam","steep riffle","gateway riffle","eye riffle","herringer riffle",
        "live oak","sunset pumps","hallwood","tisdale","knights landing")

uniqsite=siteO;Nsites=length(uniqsite)

for(isite in 1:Nsites){

  d1=subset(d0,site==uniqsite[isite])
  DoSite=uniqsite[isite]

  uniqyr=unique(sort(d1$run_year))
  Nyrs=length(uniqyr)

  if(DoSite=="ucc"){ #load pcap fit object needed by prepare_abundance_inputs. Sites ordered by trib then mainstem so only load trib pcap object for first trib
    load(file="C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_trib.rdata")
  } else if (DoSite=="tisdale" | DoSite=="knights landing"){
    load(file=paste0("C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/Output/pCap_mainstem_skew_re_",DoSite,".rdata"))
  }

  for(iyr in 1:Nyrs){#1:Nyrs
    DoYr=uniqyr[iyr]

    abundance_inputs <- prepare_abundance_inputs(site = DoSite, run_year = DoYr, effort_adjust=T,minPcapMult=0.5, pcap_model_object = pcap)

    abundance <- fit_abundance_model_BUGS(abundance_inputs,
                                          "C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/model_files/abundance_model.bug",
                                          "C:/Projects/BayDelta/SAC_JPE/SRJPEmodel/data-raw/WinBUGS14")

    save(abundance,file=paste0(OutDir,"/",DoSite,"_",DoYr,".Rdata"))

  }
}
