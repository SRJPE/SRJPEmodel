library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

source("GetData.R")

#load(file="Results/fit_NoCov.Rdata")
#load(file="Results/fit_Year_Reach.Rdata")
#load(file="Results/fit_CovWY2.Rdata")
#load(file="Results/fit_CovWY2_Reach.Rdata")
load(file="Results/fit_CovWY3_Reach.Rdata")


par(mfcol=c(1,1),mai=c(0.85,0.85,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.3)#cex.cor * r
}

#Pair plots of Fixed effects
#parlist=c("S_bReach","S_bCov","SurvWoodSac","SurvForecast","S_RE","RE_sd")#"pred_surv",
#S_bReach=as.data.frame(fit, pars = c("S_bReach"))
S_bCov=as.data.frame(fit, pars = c("S_bCov"))
RE_sd=as.data.frame(fit, pars = c("RE_sd"))
xpars=cbind(S_bReach,S_bCov,RE_sd)

#If pair plotting pcap and survival for a particular release group or year
pred_pcap=as.data.frame(fit, pars = c("pred_pcap"));pred_surv=as.data.frame(fit, pars = c("pred_surv"))
Nsims=dim(pred_surv)[1];pcap=matrix(nrow=Nsims,ncol=Nreaches);surv=pcap

irg=3
#iyr=irg #for year_reach model. In this case both pcap and survival are estimated by year so set irg to the year (1-8_)
iyr=switch(irg,1,1,2,3,4,4,4,5,5,6,7,7,8) #for an irg select the appropriate index for the year for pCap

for(j in 1:Nreaches){
    vname=paste0("pred_pcap[",iyr,",",j,"]");icol=which(names(pred_pcap)==vname);  pcap[,j]=pred_pcap[,icol]
    vname=paste0("pred_surv[",irg,",",j,"]");icol=which(names(pred_surv)==vname);  surv[,j]=pred_surv[,icol]
}
pcap=as.data.frame(pcap);colnames(pcap)=paste0("pcap",1:Nreaches)
surv=as.data.frame(surv);colnames(surv)=paste0("surv",1:Nreaches)
#xpars=cbind(pcap,surv,RE_sd)

#If pair plotting a set of pCaps or survs from one reach over years for year_reach model
#ir=3
#for(iyr in 1:Nyrs){
#  vname=paste0("pred_pcap[",iyr,",",ir,"]");icol=which(names(pred_pcap)==vname);  pcap[,iyr]=pred_pcap[,icol]
  #vname=paste0("pred_surv[",irg,",",j,"]");icol=which(names(pred_surv)==vname);  surv[,j]=pred_surv[,icol]
#}
#pcap=as.data.frame(pcap);colnames(pcap)=paste0("pcap",1:Nyrs)
#xpars=c(pcap)

pairs(x=xpars,pch=19,cex=0.2,diag.panel=panel.hist,upper.panel=panel.cor)#,main=sname

