rm(list=ls(all=TRUE))
library("R2WinBUGS")

Nmcmc=40000;Nburnin=20000;Nchains=3;Nthin=20
OutDir="Results/"
BaseMod="Rsp_Mod"

d=read.table(file="Battle.txt",header=T)
Nyrs=dim(d)[1]

obs_redd=d$spcnt
obs_passage=d$upcnt
X0=d$covar
X=(X0-mean(X0))/sd(X0)


data<-list("Nyrs","obs_redd","obs_passage","X")	
inits1<-list(logit_RsP=rep(0,Nyrs),b1_surv=0,tauRsP=1)
inits2<-list(logit_RsP=rep(0,Nyrs),b1_surv=0,tauRsP=1)
inits3<-list(logit_RsP=rep(0,Nyrs),b1_surv=0,tauRsP=1)
inits<-list(inits1,inits2,inits3)
parameters<-c("RsP","b1_surv","pred_redd","muRsP","tauRsP","sdRsP")

ModName=paste0(BaseMod,".bug")
sr.sim<-bugs(data, inits, parameters, ModName, n.chains=Nchains, n.burnin=Nburnin, n.thin=Nthin, n.iter=Nmcmc, debug=T,codaPkg=F,DIC=T,clearWD=F,bugs.directory="c:/WinBUGS14/")

fn_post=paste(BaseMod,"post.out",sep="");write.table(file=fn_post,sr.sim$sims.list,row.names=F,col.names=T)
fn_sum=paste(BaseMod,"sum.out",sep="");write.table(file=fn_sum,sr.sim$summary,append=F)
fn_dic=paste(BaseMod,"dic.out",sep="");write(file=fn_dic,c(sr.sim$pD,sr.sim$DIC),ncolumns=2)



