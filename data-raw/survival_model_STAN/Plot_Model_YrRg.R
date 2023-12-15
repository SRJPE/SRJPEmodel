library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

source("GetData.R")
#load(file="Results/fit_Year_Reach.Rdata");N=Nyrs;xticklab=Year;Nfor=0;   rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0)
load(file="Results/fit_YR_HBM.Rdata");N=Nyrs;xticklab=Year;Nfor=0;   rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0)
#load(file="Results/fit_CovWY2.Rdata");N=Nrg;xticklab=RelGp;Nfor=2;ForLab=c("Dry","Wet");   rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0)
#load(file="Results/fit_CovWY2_Reach.Rdata");N=Nrg;xticklab=RelGp;Nfor=2;ForLab=c("Dry","Wet");  rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0)
#load(file="Results/fit_CovWY3_Reach.Rdata");N=Nrg;xticklab=RelGp;Nfor=3;ForLab=c("Critical","D/BN","Wet"); rgwy_ind=c(1,1,0,1,2,2,2,1,1,2,1,1,0)

ncol=length(unique(rgwy_ind))
if(ncol==2){
  xcol=c("red","blue")
} else {
  xcol=c("red","grey","blue")
}

#par(mfcol=c(2,2),mai=c(1,0.8,0.3,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1,font.lab=1)
par(mfcol=c(2,1),mai=c(0.7,0.8,0.2,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1,font.lab=1)
#par(mfcol=c(1,1),mai=c(1.2,0.8,0.3,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1,font.lab=1)

Plot_Mat<-function(N,varnm,ylabel,xticklab,pleg){
	Nx=N*Nreaches+N
	xcol2=c("blue","green","orange","red")

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
	    stats[2]=mean(mat[,icol])
	    k=k+1
	    points(k,stats[2],pch=19,cex=1.2,col=xcol2[j])
	    arrows(x0=k,x1=k,y0=stats[1],y1=stats[3],col=xcol2[j],angle=90,code=3,length=0.025)
	    if(j==2)xtick_loc[i]=k
	  }
	  k=k+1
	  abline(v=k,lty=2,col="grey")
	}
	axis(1,at=xtick_loc,labels=xticklab,las=3,cex.axis=0.8)
	if(pleg==T) legend("bottomleft",horiz=F,legend=c("Woodson","Butte","Sacramento","Delta"),pch=rep(19,4),col=xcol2,bty='n',ncol=2)
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
		stats[i,2]=mean(survmat[,icol])
	}
	plot(1:N,stats[,2],type='p',pch=19,cex=1.2,ylim=c(0,max(stats)),bty='l',xlab="",ylab=ylabel,axes=F)
	axis(2);axis(1,at=1:N,labels=xticklab,las=3,cex.axis=0.8)
	
	for(i in 1:N) {
	  points(x=i,y=stats[i,2],col=xcol[rgwy_ind[i]+1],pch=19,cex=1.2)
	  #lines(x=c(i,i),y=c(stats[i,1],stats[i,3]),col=xcol[rgwy_ind[i]+1])
	  arrows(x0=i,x1=i,y0=stats[i,1],y1=stats[i,3],col=xcol[rgwy_ind[i]+1],angle=90,code=3,length=0.025)
	}
	legend("topleft",legend=ForLab,col=xcol,pch=rep(19,ncol),bty='n')

}

Plot_Vec(N,varnm="SurvWoodSac",ylabel="Woodson-Sacramento Survival",xticklab)

#Plot forecast
if(Nfor>0){
	varnm="SurvForecast"
	sfor=as.data.frame(fit,pars=varnm)
	stats=matrix(nrow=Nfor,ncol=3)
	for(i in 1:Nfor){
		vname=paste0(varnm,"[",i,"]")
		icol=which(names(sfor)==vname)
		#stats[i,]=quantile(sfor[,icol],prob=c(0.025,0.5,0.975))
		stats[i,]=quantile(sfor[,icol],prob=c(0.16,0.5,0.84))
		stats[i,2]=mean(sfor[,icol])
		#stats[i,1]=stats[i,2]-sd(sfor[,icol]);stats[i,3]=stats[i,2]+sd(sfor[,icol])
		
	}
	plot(1:Nfor,stats[,2],type='p',pch=19,cex=1.2,ylim=c(0,max(stats)),bty='l',xlab="",ylab="Woodson-Sacramento Survival Forecast",axes=F)
	axis(2);axis(1,at=1:Nfor,labels=ForLab,las=3,cex.lab=1)
	for(i in 1:Nfor){
	  #lines(x=c(i,i),y=c(stats[i,1],stats[i,3]))
	  points(i,stats[i,2],pch=19,cex=1.2,col=xcol[i])
	  arrows(x0=i,x1=i,y0=stats[i,1],y1=stats[i,3],col=xcol[i],angle=90,code=3,length=0.025)
	}
}

