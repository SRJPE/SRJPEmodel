#Cycle through a series of tribs and years hat have special priors and run model
rm(list=ls(all=TRUE))

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

MultiRun_Mode=T
doWks=c(seq(45,53),seq(1,22))

#Identify unique combinations of trib_rst sites and run years in special prior file
dsp=read.csv(file="data-raw/juvenile_abundance/btspas_model_code/Special_Priors.csv",header=T)
uniq_tr=unique(dsp[c("Stream_Site","RunYr")])
Nruns=dim(uniq_tr)[1]

for(ii in 1:Nruns){
	doTrib=uniq_tr[ii,1]
	doYr=uniq_tr[ii,2]

	source("data-raw/juvenile_abundance/btspas_model_code/Model.R")
	#source("PlotModel.R")

}
