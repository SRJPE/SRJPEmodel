
#graphics.off()
#if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
#windows(record=T,width=12,height=12);par(las=1)

fnm=paste0("Output/Data.Rdata")
load(file=fnm)#load model_data list

fnm="Output/Fit.Rdata"#read in fit from previous run
load(file=fnm)

Vnm=c("pred_lgNtot","pred_lgNx","lt_cp","JPE_trib")
dp=as.data.frame(fit, pars = Vnm)#;Nsims=dim(dp)[1]

  
  for(itrib in 1:Ntribs){
    par(mfcol=c(3,4),xaxs="i",mai=c(.8,.7,.1,.1), omi=c(0.25,0.25,0.5,0.25))
    
    for(ifor in 1:Nfor){
    
      for(ivar in 1:3){
        
        if(ivar==1){
          Vlab="Log Outmigrant Abundance at RST ('000s)"
          Vline=-50
        } else if (ivar==2){
          Vlab="Log Outmigrant Abundance through Forecast Week ('000s)"
          Vline=-27
        } else if (ivar==3){
          Vlab="Logit Proportion of Outmigrants Past RST through Forecast Week"
          Vline=-2
        }
        
        if(ivar==1){
          prmu=model_data$srNtot_mu[itrib];prsd=model_data$srNtot_sd[itrib]  #doesn't change across forecast wks but must be in loop because pr values for other ivar's do
        } else if (ivar==2){
          prmu=model_data$obs_Nx_mu[itrib,ifor];prsd=model_data$obs_Nx_sd[itrib,ifor] 
        } else if (ivar==3){
          prmu=model_data$cp_mu[itrib,ifor];prsd=model_data$cp_sd[itrib,ifor] 
        }
        
        icol=which(names(dp)==paste0(Vnm[ivar],"[",itrib,",",ifor,"]"))
        post=dp[,icol] 
        
        
        prior=rnorm(n=10000,mean=prmu,sd=prsd)
        Pr=density(prior)
        
        Po=hist(post,plot=F)  
        ymax=max(c(max(Po$density),max(Pr$y)))
        
        xmain=For_CalWk[ifor]
        hist(post,ylim=c(0,ymax),xlab="",ylab="",main=xmain,freq=F,axes=F)
        axis(1)
        lines(density(prior),lty=2)
        if(ivar==1 & ifor==1) legend("topright",legend=c("posterior","prior"),lty=c(NA,2),pch=c(22,NA),pt.bg=c("gray","NA"),bty='n')
      
        if(ifor==1){
          #text(x=min(Po$mids),y=ymax*0.8,adj=0,labels=Vlab,cex=0.5,font=2)
          mtext(Vlab,side=1, las=1,line=Vline,outer=T,cex=1.1,font=2)
        }
      }#next ivar    
    }#next ifor
      
    
    mtext("Frequency",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)
    mtext(DoTrib[itrib],side=3, las=1,line=1,outer=T,cex=1.3,font=2)
  
    
}#next itrib

