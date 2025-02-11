#Cycle through a series of tribs and years and run model
rm(list=ls(all=TRUE)) 

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

library("R2WinBUGS")
library(splines2)

MultiRun_Mode=T

#pdf(file="AllPlots.pdf")	#if writing plots to a file

dr0a=read.table(file="RST_Input.txt",header=T)
dr0=subset(dr0a,Trib!="sacramento river_knights landing" & Trib!="sacramento river_tisdale")

dr=unique(dr0[c("Trib","RunYr")])	#list of unique tribs and run years
Nruns=dim(dr)[1]			# # of runs

doWks=c(seq(45,53),seq(1,22)) #can't have a wider range than allwks

MinSampProp=0.75	

Nests=0
	 
for(ii in 1:Nruns){
	doTrib=dr$Trib[ii]
	doYr=dr$RunYr[ii]
	
	dsc=subset(dr0,Trib==doTrib & RunYr==doYr & is.na(ua)==F)
	nwks=dim(dsc)[1]
	
	if(nwks/length(doWks)>=MinSampProp){ #Don't process if number of weeks for Trib-RunYr < min sample propprtion
		#source("Model.R")
		Nests=Nests+1
		source("PlotModel.R")
	}
}
#dev.off()#if writing plots to a file.