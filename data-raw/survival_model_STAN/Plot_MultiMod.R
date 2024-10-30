#Plot results from multiple models on the same page

library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

source("GetData.R")
yrline=c(2.5,3.5,4.5,7.5,9.5,10.5,12.5)

Plot_Vec<-function(N,varnm,ymax,xticklab,xmain){
  
  survmat=as.data.frame(fit, pars = c(varnm))
  stats=matrix(nrow=N,ncol=3)
  for(i in 1:N){
    vname=paste0(varnm,"[",i,"]")
    icol=which(names(survmat)==vname)
    stats[i,]=quantile(survmat[,icol],prob=c(0.025,0.5,0.975))
    stats[i,2]=mean(survmat[,icol])
  }
  plot(1:N,stats[,2],type='p',pch=19,cex=1.2,ylim=c(0,ymax),bty='l',xlab="",ylab="",axes=F,main=xmain)
  axis(2);axis(1,at=1:N,labels=xticklab,las=3,cex.axis=0.6)
  
  for(i in 1:N) {
    points(x=i,y=stats[i,2],col=xcol[rgwy_ind[i]+1],pch=19,cex=1.2)
    #lines(x=c(i,i),y=c(stats[i,1],stats[i,3]),col=xcol[rgwy_ind[i]+1])
    arrows(x0=i,x1=i,y0=stats[i,1],y1=stats[i,3],col=xcol[rgwy_ind[i]+1],angle=90,code=3,length=0.025)
  }
  abline(v=yrline,col="grey",lty=2)
  legend("topleft",legend=ForLab,col=xcol,pch=rep(19,ncol),bty='n')
  
}

Plot_For<-function(Nfor,varnm,ymax,ForLab,xmain){

  sfor=as.data.frame(fit,pars=varnm)
  stats=matrix(nrow=Nfor,ncol=3)
  for(i in 1:Nfor){
    if(Nfor>1){
        vname=paste0(varnm,"[",i,"]")
    } else {
        vname=varnm
    }
    
    icol=which(names(sfor)==vname)
    #stats[i,]=quantile(sfor[,icol],prob=c(0.025,0.5,0.975))
    stats[i,]=quantile(sfor[,icol],prob=c(0.16,0.5,0.84)) #1SD which is 16-84th percentiles
    #stats[i,]=quantile(sfor[,icol],prob=c(0.25,0.5,0.75))
    stats[i,2]=mean(sfor[,icol])
  }

  plot(1:Nfor,stats[,2],type='p',pch=19,cex=1.2,ylim=c(0,ymax),bty='l',xlab="",ylab="",main=xmain,axes=F)
  axis(2);axis(1,at=1:Nfor,labels=ForLab,las=3,cex.lab=1)
  for(i in 1:Nfor){
    #lines(x=c(i,i),y=c(stats[i,1],stats[i,3]))
    points(i,stats[i,2],pch=19,cex=1.2,col=xcol[i])
    arrows(x0=i,x1=i,y0=stats[i,1],y1=stats[i,3],col=xcol[i],angle=90,code=3,length=0.025)
  }
  abline(h=seq(0.1,ymax-0.1,by=0.1),col="grey",lty=2)
}


par(mfrow=c(2,2),mai=c(.7,0.8,0.1,0.1), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1,font.lab=1)

for(imod in 1:4){
  if(imod==1){
    load(file="Results/fit_NoCov.Rdata");xmain="Null";N=Nrg;xticklab=rep("",Nrg);Nfor=1;ForLab="All";   rgwy_ind=rep(0,Nrg)
  } else if (imod==2){
    load(file="Results/fit_CovWY2.Rdata");xmain="2 WY Types";N=Nrg;xticklab=rep("",Nrg);Nfor=2;ForLab=c("Dry","Wet");   rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0)
  } else if (imod==3){
    load(file="Results/fit_CovWY2_Reach.Rdata");xmain="2 WY Types*Reach";N=Nrg;xticklab=RelGp;Nfor=2;ForLab=c("Dry","Wet");  rgwy_ind=c(0,0,0,1,1,1,1,1,1,1,0,0,0)
  } else if (imod==4){
    load(file="Results/fit_CovWY3_Reach.Rdata");xmain="3 WY Types*Reach";N=Nrg;xticklab=RelGp;Nfor=3;ForLab=c("Critical","D/BN/AN","Wet"); rgwy_ind=c(1,1,0,1,2,2,2,1,1,2,1,1,0)
  }
  
  ncol=length(unique(rgwy_ind))
  if(ncol==1){
    xcol="black"
  } else if(ncol==2){
    xcol=c("red","blue")
  } else {
    xcol=c("red","grey","blue")
  }

  Plot_Vec(N,varnm="SurvWoodSac",ymax=0.6,xticklab,xmain)
  #Plot_For(Nfor,varnm="SurvForecast",ymax=0.6,ForLab,xmain)
}
#mtext("Water Year Type",side=1, outer=T,,cex=1.2,font=2,line=1)
mtext("Woodson-Sacramento Survival Rate",side=2, las=3,outer=T,,cex=1.2,font=2,line=-1)