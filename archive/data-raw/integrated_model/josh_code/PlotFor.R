#Plot outmigrants by trib and in total based on SR, inseason, and stan model. 
par(mfrow=c(2,3),xaxs="i",mai=c(.8,.7,.5,.5), omi=c(0.25,0.25,0.5,0.25),cex.axis=0.9)
#par(mfrow=c(1,1),xaxs="i",mai=c(1,1,.1,.1), omi=c(0.5,0.5,0.5,0.5))
#for(itrib in 1:1){
for(itrib in 1:Ntribs){
  #ymax=max(predSR_stats[itrib,3],predIn_stats[itrib,,3],pred_Ntot_stats[itrib,ifor,3])
  ymax=max(predSR_stats[itrib,3],predIn_stats[itrib,3:Nfor,3])
  #if(ymax>2*predSR_stats[itrib,3]) ymax=max(predIn_stats[itrib,,2])*1.5
  x=barplot(height=c(predSR_stats[itrib,2],predIn_stats[itrib,,2]),las=3,names.arg=c("SR",For_CalWk),main=DoTrib[itrib],ylim=c(0,ymax))
  i=1; arrows(x0=x[i],x1=x[i],y0=predSR_stats[itrib,1],y1=predSR_stats[itrib,3],angle=90,code=3,length=0.025,col="black")
  for(ifor in 1:Nfor){
    arrows(x0=x[ifor+1],x1=x[ifor+1],y0=predIn_stats[itrib,ifor,1],y1=predIn_stats[itrib,ifor,3],angle=90,code=3,length=0.025,col="black")
    points(x=x[ifor+1]+0.1,y=pred_Ntot_stats[itrib,ifor,2],pch=19,col="red")
    arrows(x0=x[ifor+1]+0.1,x1=x[ifor+1]+0.1,y0=pred_Ntot_stats[itrib,ifor,1],y1=pred_Ntot_stats[itrib,ifor,3],angle=90,code=3,length=0.025,col="red")
    
  }
  if(itrib==1){
    legend("topright",legend=c("Independent","Joint SR-Inseason"),pch=c(15,19),pt.cex=rep(1.5,2),col=c("grey","red"),bty='n')
  }
}
mtext("Outmigrant abundance ('000s)",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)


#summed across all tribs
par(mfrow=c(1,1),xaxs="i",mai=c(.8,.8,.7,.7), omi=c(0.25,0.25,0.5,0.25),cex.lab=1.3,cex.main=1.3,cex.axis=1.3)
xoff=0.075
ymax=max(pred_Ntot_all_stats)
plot(x=1:Nfor-xoff,y=pred_Ntot_all_stats[,2],type='p',bty='n',pch=19,cex=1.5,main="All Tributaries Combined",ylim=c(0,ymax),xlim=c(0.5,Nfor+0.5),col="red",xlab="",ylab="",axes=F)
axis(2);axis(1,at=1:Nfor,labels=For_CalWk,las=3,cex=1.2)
mtext("at Rotary Screw Traps (000s)", side=2, line=4, cex=1.3,las=3, font=2,col="red")
for(ifor in 1:Nfor) arrows(x0=ifor-xoff,x1=ifor-xoff,y0=pred_Ntot_all_stats[ifor,1],y1=pred_Ntot_all_stats[ifor,3],angle=90,code=3,length=0.025,col="red")

par(new=T)
plot(x=1:Nfor+xoff,y=JPE_stats[1:Nfor,2],,xlim=c(0.5,Nfor+0.5),ylim=c(0,ymax*0.1),type='p',bty='n',xlab="",ylab="",pch=19,cex=1.5,col="blue",axes=F)
axis(4);mtext("at Delta Entry ('000s)", side=4, line=3, cex=1.3,las=3, font=2,col="blue")
for(ifor in 1:Nfor) arrows(x0=ifor+xoff,x1=ifor+xoff,y0=JPE_stats[ifor,1],y1=JPE_stats[ifor,3],angle=90,code=3,length=0.025,col="blue")
#legend("topleft",legend=c("at RST","at Delta"),pch=c(19,19),pt.cex=rep(1.5,2),col=c("red","blue"),bty='n')


#Downstream survival rate
par(mfrow=c(1,1),xaxs="i",mai=c(1,1,.7,.7), omi=c(0.25,0.25,0.5,0.25),cex.lab=1.3,cex.main=1.3,cex.axis=1.3)
ymax=max(pred_DSsurv_stats)
x=barplot(pred_DSsurv_stats[,2],ylim=c(0,ymax),names.arg=DoTrib,bty='n',xlab="",ylab="")
for(itrib in 1:Ntribs) arrows(x0=x[itrib],x1=x[itrib],y0=pred_DSsurv_stats[itrib,1],y1=pred_DSsurv_stats[itrib,3],angle=90,code=3,length=0.025)
mtext("RST-Delta Entry Survival Rate", side=2, line=4, cex=1.3,las=3, font=2)
mtext("RST Site", side=1, line=3, cex=1.3, font=2)
