#rm(list=ls(all=TRUE))

library(scales)
lb=0.025;ub=0.975	#0.975
Nscale=0.001		#abundance is in units of '000s of fish

dx=read.table(file="../Data/Jwk_Dates.txt",header=T)
Jdt=dx$Date

PlotAbOnly=T
PlotIt=T

MultiRun_Mode=T
if(MultiRun_Mode==F){
	doTrib="battle creek_ubc";doYr=2009
	#doTrib="deer creek_deer creek";doYr=2010
	#doTrib="clear creek_lcc";doYr=2005
	#doTrib="feather river_eye riffle";doYr=2020
	#doTrib="butte creek_okie dam 1";doYr=2007
	
	#RunNm=paste0("output/",doTrib,"_",doYr)
	RunNm=paste0(doTrib,"_",doYr)
	#RunNm=""
	outdir=""
	
	graphics.off()
	if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
	windows(record=T);par(las=1)
	
} else{
	#LastPos=as.integer(regexpr(pattern=" ",doTrib))-1;RunNm=paste0(strtrim(doTrib,width=LastPos),"_",doYr)
	RunNm=paste0(doTrib,"_",doYr)
	outdir="output/"
	outdir2="OutSpecPriors/"
	
	#check to see if there is a matching file in OutSpecPriors and if so, swap the directory so it is read in and plotted
	fn0=paste0(RunNm,"_sum.out")
	file.ls2=intersect(list.files(path=outdir2),list.files(path=outdir2,pattern="sum.out"))
	if(is.na(match(fn0,file.ls2))==F) outdir=outdir2
	
	#check to see if the trib-year was modelled and therefore has *sum.out file. If it doesn't exist don't plot
	file.ls1=list.files(path=outdir,pattern="sum.out")
	if(is.na(match(fn0,file.ls1))==T) PlotIt=F
}

if(PlotIt==T){
	fn_data=paste0(outdir,RunNm,"_data.out")
	d=read.table(file=fn_data,header=T)
	Nstrata=dim(d)[1]

	fn_pred=paste0(outdir,RunNm,"_post.out")
	dp=read.table(file=fn_pred,header=T)
	Nsims=dim(dp)[1]

	fn_sum=paste0(outdir,RunNm,"_sum.out")
	dsum=read.table(file=fn_sum,header=T)
	irow=which(rownames(dsum)=="Ntot")
	Ntot_Rhat=dsum$Rhat[irow]

	fn_knots=paste0(outdir,RunNm,"_knots.out")
	kknots=scan(file=fn_knots)
	Nknots=length(kknots)


	dmr=read.table(file="MRdata.txt",header=T)
	Ntribs=length(unique(dmr$Trib))
	Trib=character(length=Ntribs)
	for(i in 1:Ntribs)Trib[i]=dmr$Trib[min(which(dmr$TribInd==i))]

	if(PlotAbOnly==F){
		par(mfcol=c(2,1),mai=c(0.8,1,0.2,0.7), omi=c(0.1,0.1,0.3,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
	} else {
		par(mfcol=c(1,1),mai=c(1,1,0.2,0.7), omi=c(0.1,0.1,0.3,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)
	}

	####Plot weekly abundance ####################
	Ntot=as.double(quantile(dp$Ntot,prob=c(lb,0.5,ub)))*Nscale
	Ntot_cv=round(100*sd(dp$Ntot)/mean(dp$Ntot),digits=0)
	NtotLab=paste0(round(Ntot[2],digits=0)," (",round(Ntot[1],digits=0)," - ", round(Ntot[3],digits=0), ") cv=",Ntot_cv,"%")

	N=matrix(nrow=Nstrata,ncol=3)
	for(i in 1:Nstrata){
		vn=paste0("N.",i);icol=which(names(dp)==vn)
		N[i,]=as.double(quantile(dp[,icol]*Nscale,prob=c(lb,0.5,ub)))
	}

	#Lincoln_Peterson predictions of abundance by strata
	irecs_N=which(is.na(d$Releases)==F & d$Recaptures>0 & is.na(d$u)==F)
	pN=(d$u[irecs_N]/(d$Recaptures[irecs_N]/d$Releases[irecs_N]))*Nscale
	sdN=sqrt(Nscale^2*(as.numeric(d$Releases[irecs_N])*as.numeric(d$Recaptures[irecs_N])*(as.numeric(d$Releases[irecs_N])-as.numeric(d$Recaptures[irecs_N]))*d$u[irecs_N])/(d$Recaptures[irecs_N]^2*d$Recaptures[irecs_N]))

	#Chapman estimate which is LP+1 stuff
	#pN=(((d$Releases[irecs_N]+1)*(d$u[irecs_N]+1))/(d$Recaptures[irecs_N])-1)*Nscale
	#sdN=sqrt(Nscale^2*((d$Releases[irecs_N]+1)*(d$Recaptures[irecs_N]+1)*(d$Releases[irecs_N]-d$Recaptures[irecs_N])*d$u[irecs_N])/((d$Recaptures[irecs_N]+1)^2*(d$Recaptures[irecs_N]+2)))

	xmain=paste0(RunNm," Ntot=",NtotLab)
	ylims=c(-0.1,max(c(N,pN))*1.1)
	irecs=which(is.na(d$u)==T);#points(irecs,rep(ylims[1],length(irecs)),pch=19,cex=1.1,col='red') #identify strata that were not sampled
	bcols=rep("grey",Nstrata);bcols[irecs]="red"
	x=barplot(height=N[,2],col=bcols,names.arg=Jdt[d$Jwk],cex.names=0.7, las=3,axis.lty=1,ylim=ylims,  bty='l',main=xmain, xlab="",ylab="Abundance ('000s)")
	#points(x=x[round(kknots,digits=0)],y=rep(ylims[2]*0.8,Nknots),pch=25,bg="green",cex=1.5)#knot position

	for(i in 1:Nstrata){#error bars and catch
		arrows(x0=x[i],x1=x[i],y0=N[i,1],y1=N[i,3],angle=90,code=3,length=0.025,col="black")
		par(srt=90);text(x=x[i],y=ylims[2]*0.9,labels=d$u[i],cex=0.7);par(srt=0)
		if(is.na(match(i,irecs))==F) points(x[i],ylims[2]*0.9,pch=19,col="red",cex=0.8)#identify strata with no sampling (red bars may not show up if estimated abundance is very low)
	}

	#Plot Peterson estimates and CL's
	coff=0.1
	nn=length(irecs_N)
	if(nn>0){
		for(i in 1:nn){
			points(x[irecs_N[i]]-coff,pN[i],pch=21,cex=1.5,col='blue')	#data-driven predictions of N
			if(sdN[i]>0){
				lc=pN[i]-1.96*sdN[i];uc=pN[i]+1.96*sdN[i]
				arrows(x0=x[irecs_N[i]]-coff,x1=x[irecs_N[i]]-coff,y0=lc,y1=uc,angle=90,code=3,length=0.025,col="blue")
			}
		}
	}
	#legend("left",legend=c("Bayesian","Peterson","No sampling","Spline knot"),pch=c(22,21,19,25),pt.cex=rep(1.2,4),pt.bg=c("grey",NA,NA,"green"),col=c("grey","blue","red","green"),cex=0.8,bty='n')
	legend("right",legend=c("Bayesian","Peterson","No sampling or no releases","Discharge"),pch=c(22,21,19,NA),pt.cex=c(rep(1.2,3),NA),pt.bg=c("grey",NA,NA,NA),col=c("grey","blue","red","darkgreen"),lwd=c(NA,NA,NA,1.2),cex=0.8,bty='n')
	
	if(PlotAbOnly==T){
		par(new=T)
		plot(x-coff,d$Flow/1000,ylim=c(min(d$Flow/1000),max(d$Flow/1000)*1.2), xlim=c(0,trunc(max(x))+1),type='o',pch=19,cex=0.7,lwd=2,col="darkgreen",axes=F,xlab="",ylab="");axis(4,cex.axis=0.9)
		mtext("Discharge ('000s cfs)",side=4, las=3,adj=0.5,outer=T,cex=1.0,font=2,line=-1,col="darkgreen")
	}
	
	if(PlotAbOnly==F){
		###Plot weekly capture probability #########################
		pCap=matrix(nrow=Nstrata,ncol=3)
		for(i in 1:Nstrata){
			vn=paste0("pCap_U.",i);icol=which(names(dp)==vn)
			pCap[i,]=as.double(quantile(dp[,icol],prob=c(lb,0.5,ub)))
		}
		ylims=c(0,max(pCap)*1.1)
		irecs=which(is.na(d$Releases)==T)#;points(irecs,pCap[irecs,2],pch=19,cex=1.1,col='red')
		bcols=rep("grey",Nstrata);bcols[irecs]="red"
		x=barplot(height=pCap[,2],col=bcols,names.arg=Jdt[d$Jwk],cex.names=0.7, las=3,axis.lty=1,ylim=ylims,  bty='l', xlab="First Date of Week",ylab="Capture Probability")

		irecs=which(is.na(d$Releases)==F)
		points(x[irecs]-coff,d$Recaptures[irecs]/d$Releases[irecs],pch=21,cex=1.5,col='blue')

		for(i in 1:Nstrata) {
			arrows(x0=x[i],x1=x[i],y0=pCap[i,1],y1=pCap[i,3],angle=90,code=3,length=0.025,col="black")
			par(srt=90)
			text(x=x[i],y=ylims[2]*0.95,labels=d$Recaptures[i],cex=0.7)
			text(x=x[i],y=ylims[2]*0.82,labels=d$Releases[i],cex=0.7)
			par(srt=0)
		}

		par(new=T)
		plot(x-coff,d$Flow/1000,ylim=c(min(d$Flow/1000),max(d$Flow/1000)*1.2), xlim=c(0,trunc(max(x))+1),type='o',pch=19,cex=0.7,lwd=2,col="darkgreen",axes=F,xlab="",ylab="");axis(4,cex.axis=0.9)
		mtext("                           Discharge ('000s cfs)",side=4, las=3,adj=0,outer=T,,cex=1.0,font=2,line=-1,col="darkgreen")
	}
		print(c(doTrib,doYr,Ntot_Rhat))

		### Plot pCap distributions ####
		if(MultiRun_Mode==F){

			#pCap distributions for each tributary in one panel and overlay hyper distribution
			par(mfcol=c(1,1),mai=c(0.85,2.25,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
			xmax=0.3
			pstat=matrix(nrow=Ntribs,ncol=3)
			for(i in 1:Ntribs){
				vn=paste0("b0_pCap.",i);icol=which(names(dp)==vn)
				pCap=exp(dp[,icol])/(1+exp(dp[,icol]))
				pstat[i,]=as.double(quantile(pCap,prob=c(0.025,0.5,0.975)))
			}
			plot(pstat[,2],1:Ntribs,pch=19,bty='l',xlim=c(0,xmax),axes=F,ylab="",xlab="Mean Trap Efficiency")

			tnm=character(length=Ntribs);mup=vector(length=Ntribs)
			xTrib=Trib#make a copy since subset gets confused when dataframe variable name is the same as the one it is checking against
			mup=vector(length=Ntribs)
			for(i in 1:Ntribs){
				dmr1=subset(dmr,Trib==xTrib[i] & is.na(Rel1)==F)
				mup[i]=sum(dmr1$Recap1)/sum(dmr1$Rel1)	#mup[i]=mean(dmr1$Recap1/dmr1$Rel1)
				nsamps=dim(dmr1)[1]
				tnm[i]=paste0(gsub(pattern="riffle",replacement="",x=Trib[i])," (",nsamps,")")
			}
			axis(1);axis(2,at=1:Ntribs,labels=tnm,cex.axis=0.75)

			for(i in 1:Ntribs){#plot CI's and point estimate for each trib
				arrows(x0=pstat[i,1],x1=pstat[i,3],y0=i,y1=i,angle=90,code=3,length=0.025,col="black")
				points(x=mup[i],y=i,pch=21,cex=1.5,col='blue')
			}
			hmu=mean(dp$trib_mu.P);	hmu=exp(hmu)/(1+exp(hmu)); abline(v=hmu,lty=2)

			#hyper distribution
			pvec=seq(0,xmax,length.out=50);logit_pvec=log(pvec/(1-pvec))
			prob=dnorm(x=logit_pvec,mean=mean(dp$trib_mu.P),sd=mean(dp$trib_sd.P))
			par(new=T); plot(pvec,prob,lty=1,type='l',xlim=c(0,xmax),axes=F,xlab="",ylab="")
			legend("topright",pch=c(19,21,NA),lty=c(NA,NA,1),col=c("black","blue","black"),legend=c("Model Estimate","Point Estimate","Hyper-distribution"),bty='n')


			par(mfcol=c(1,2),mai=c(1,0.75,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
			#hyper-distributions from which pCaps will be drawn for tribs with no efficiency estimates
			b0=vector(length=Nsims);b1=b0
			for(i in 1:Nsims){
				b0[i]=rnorm(n=1,dp$trib_mu.P[i],dp$trib_sd.P[i])
				b1[i]=rnorm(n=1,dp$flow_mu.P[i],dp$flow_sd.P[i])
			}
			pCap3=exp(b0)/(1+exp(b0))
			hist(pCap3,col="grey",xlab="Trap Efficiency for Tribs Without MR at Average Flow",breaks=seq(0,1,by=0.005),xlim=c(0,xmax),main="",yaxt='n',ylab="")


			#Plot discharge-trap efficiency relationship for a trib without MR data
			Q=seq(-2,2,by=0.01)
			pstat=matrix(nrow=length(Q),ncol=3)
			for(i in 1:length(Q)){
				pCap4=exp(b0+b1*Q[i])/(1+exp(b0+b1*Q[i]))
				pstat[i,]=as.double(quantile(pCap4,probs=c(lb,0.5,ub)))
				pstat[i,2]=mean(pCap4)
			}

			plot(Q,pstat[,2],type='l',lty=1,bty='l',ylim=c(0,max(pstat)),ylab="Trap Efficiency",xlab="Standardized Discharge")
			polygon(x = c(rev(Q), Q), y = c(rev(pstat[,1]), pstat[,3]), col = alpha("grey", 0.3), border = NA)

	}

}