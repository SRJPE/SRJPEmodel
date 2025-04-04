library(SRJPEdata)
rm(list=ls(all=TRUE))

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

OutDir="Output/"
JuvDir="../RST/RunSize/StanVersion/Output_Trib_SR/"

dA0=observed_adult_input #all adult data from this table in SRJPEdata

for(itype in 1:8){
  
  DoSite=switch(itype,"ubc","ubc","lcc","lcc","mill creek","deer creek","okie dam","okie dam")
  if(DoSite=="ubc"){
    DoStream="battle creek"
  } else if (DoSite=="lcc"){
    DoStream="clear creek"
  } else if (DoSite=="mill creek"){
    DoStream="mill creek"
  } else if (DoSite=="deer creek"){
    DoStream="deer creek"
  } else if (DoSite=="okie dam"){
    DoStream="butte creek"
  }
  
  Am=switch(itype,"rd","us","rd","us","rd","hd","cs","hd")
  if(Am=="rd"){
    Amethod="redd_count"
  } else if (Am=="us"){
    Amethod="upstream_estimate"
  } else if (Am=="hd"){
    Amethod="holding_count"
  } else if (Am=="cs"){
    Amethod="carcass_estimate"
  }

  #Make SR table. Begin by getting list of years with outmigrant abundance estimates
  file.ls=intersect(list.files(path=JuvDir,pattern=paste0("N_",DoSite)),list.files(path=JuvDir,pattern=".rdata"))
  NJuvYrs=length(file.ls)
  
  dA=subset(dA0,stream==DoStream & data_type==Amethod)
  NAdYrs=dim(dA)[1]
  print(c(DoSite,DoStream, Am,NJuvYrs,NAdYrs))
  #if(DoStream=="deer creek") browser()
}

#some extra info that doesn't fit nicely into mdoelleind loop
DoSite="knights landing"
JuvDir="../RST/RunSize/StanVersion/Output_Main_SR/"
file.ls=intersect(list.files(path=JuvDir,pattern=paste0("N_",DoSite)),list.files(path=JuvDir,pattern=".rdata"))
NJuvYrs=length(file.ls)
