library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

RunMainstem=T
OutDir="Output/NoCovEffect/"
#OutDir="Output/WithCovEffect/"
#OutDir="Output/lg1NoCovEffect/"
#OutDir="Output/lg1WithCovEffect/"

prior=rgamma(n=5000,shape=20,rate=10)

#par(mfrow=c(2,2),xaxs="i",mai=c(.6,.7,.1,.1), omi=c(0.25,0.25,0.5,0.25))
par(mfrow=c(2,1),mai=c(1,1,0.4,0.6), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
for(itr in 1:2){
  if(RunMainstem==F){
    doTrib=switch(itr,"ubc","lcc")
  } else {
    doTrib=switch(itr,"knights landing","tisdale")
  }
  
  source ("GetData.R")
  fnm=paste0(OutDir,doTrib,".Rdata")
  load(file=fnm)
  
  dp=as.data.frame(fit, pars = c("sd_pro"))
  post=dp[,1]
  
  Po=hist(post,plot=F);Pr=density(prior)
  ymax=max(c(max(Po$density),max(Pr$y)))
    
   hist(post,ylim=c(0,ymax),xlab="",ylab="",main=doTrib,freq=F)
   lines(density(prior),lty=2)
   if(itr==1) legend("topright",legend=c("posterior","prior"),lty=c(NA,2),pch=c(22,NA),pt.bg=c("gray","NA"),bty='n')
 }
 
  
