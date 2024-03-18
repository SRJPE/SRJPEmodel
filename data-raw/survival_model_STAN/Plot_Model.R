library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

source("GetData.R")
#load(file="Results/fit_Year_Reach.Rdata");N=Nyrs;xticklab=Year;Nfor=0
#load(file="Results/fit_CovWY2.Rdata");N=Nrg;xticklab=RelGp;Nfor=2;ForLab=c("Dry","Wet")
#load(file="Results/fit_CovWY2_Reach.Rdata");N=Nrg;xticklab=RelGp;Nfor=2;ForLab=c("Dry","Wet")
load(file="Results/fit_CovWY3_Reach.Rdata");N=Nrg;xticklab=RelGp;Nfor=3;ForLab=c("Critical","D/BN","Wet")

par(mfcol=c(2,2),mai=c(1,0.8,0.3,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1,font.lab=1)

Plot_Mat<-function(N,varnm,ylabel,xticklab,pleg){
	Nx=N*Nreaches+N
	xcol=c("blue","green","orange","red")

	plot(1:Nx,rep(0,Nx),type='p',cex=0,xlim=c(0,Nx),ylim=c(0,1),bty='l',xlab="",ylab=ylabel,axes=F)
	axis(2)

	mat=as.data.frame(fit, pars = c(varnm))
	xtick_loc=vector(length=Nyrs)
	k=0
	for(i in 1:N){
	  for(j in 1:Nreaches){
	    vname=paste0(varnm,"[",i,",",j,"]")
	    icol=which(names(mat)==vname)
	    stats=quantile(mat[,icol],prob=c(0.025,0.5,0.975))
	    k=k+1
	    points(k,stats[2],pch=19,cex=1.2,col=xcol[j])
	    lines(x=c(k,k),y=c(stats[1],stats[3]),lty=1,col=xcol[j])
	    if(j==2)xtick_loc[i]=k
	  }
	  k=k+1
	  abline(v=k,lty=2,col="grey")
	}
	axis(1,at=xtick_loc,labels=xticklab,las=3,cex.axis=0.8)
	if(pleg==T) legend("bottomleft",horiz=F,legend=c("Woodson","Butte","Sacramento","Delta"),pch=rep(19,4),col=xcol,bty='n',ncol=2)
}

Plot_Mat(Nyrs,varnm="pred_pcap",ylabel="Detection Probability",xticklab=Year,pleg=T)
Plot_Mat(N,varnm="pred_surv",ylabel="Survival Rate",xticklab,pleg=F)

Plot_Vec<-function(N,varnm,ylabel,xticklab){

	survmat=as.data.frame(fit, pars = c(varnm))
	stats=matrix(nrow=N,ncol=3)
	for(i in 1:N){
		vname=paste0(varnm,"[",i,"]")
		icol=which(names(survmat)==vname)
		stats[i,]=quantile(survmat[,icol],prob=c(0.025,0.5,0.975))
	}
	plot(1:N,stats[,2],type='p',pch=19,cex=1.2,ylim=c(0,max(stats)),bty='l',xlab="",ylab=ylabel,axes=F)
	axis(2);axis(1,at=1:N,labels=xticklab,las=3,cex.axis=0.8)
	for(i in 1:N) lines(x=c(i,i),y=c(stats[i,1],stats[i,3]))

}

#Plot forecast
if(Nfor>0){
	varnm="SurvForecast"
	sfor=as.data.frame(fit,pars=varnm)
	stats=matrix(nrow=Nfor,ncol=3)
	for(i in 1:Nfor){
		vname=paste0(varnm,"[",i,"]")
		icol=which(names(sfor)==vname)
		stats[i,]=quantile(sfor[,icol],prob=c(0.025,0.5,0.975))
		}
		plot(1:Nfor,stats[,2],type='p',pch=19,cex=1.2,ylim=c(0,max(stats)),bty='l',xlab="",ylab="Woodson-Sac Survival Forecast",axes=F)
		axis(2);axis(1,at=1:Nfor,labels=ForLab,las=3,cex.lab=1)
		for(i in 1:Nfor) lines(x=c(i,i),y=c(stats[i,1],stats[i,3]))

}

Plot_Vec(N,varnm="SurvWoodSac",ylabel="Woodson-Sac Survival",xticklab)
