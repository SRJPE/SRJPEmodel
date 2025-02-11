#Plot annual time series of abundance estimates

rm(list=ls(all=TRUE)) 

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)

OutDir="output/"
OutDir2="OutSpecPriors/"

doTrib=c("battle creek_ubc","clear creek_lcc","clear creek_ucc","mill creek_mill creek","deer creek_deer creek","butte creek_okie dam 1", "feather river_eye riffle", "feather river_gateway riffle","feather river_herringer riffle","feather river_steep riffle", "feather river_sunset pumps","yuba river_yuba river")#,"sacramento river_knights landing"
Ntribs=length(doTrib)

ShowMissYrs=T
ShowSpec=T
yscale=0.001

par(mfcol=c(2,2),mai=c(0.8,0.8,0.1,0.1), omi=c(0.1,0.1,0.3,0.1),cex.main=1.1,cex.lab=1.1,font.lab=2)

for(itrib in 1:Ntribs){
	file.ls=intersect(list.files(path=OutDir,pattern=doTrib[itrib]),list.files(path=OutDir,pattern="sum.out"))
	file.ls2=intersect(list.files(path=OutDir2,pattern=doTrib[itrib]),list.files(path=OutDir2,pattern="sum.out"))
	
	#This version just plots available years with no blanks for years that were not estimated
	Nyrs=length(file.ls)
	yr=character(length=Nyrs)
	mu=vector(length=Nyrs);mu=NA;mu2=mu
	cl=matrix(nrow=Nyrs,ncol=2,data=NA);cl2=cl
		
	for(i in 1:Nyrs){
		d=read.table(file=paste0(OutDir,file.ls[i]),header=T)
		irow=which(rownames(d)=="Ntot")
		mu[i]=d$mean[irow]*yscale
		cl[i,]=c(d$X2.5.[irow],d$X97.5.[irow])*yscale
		
		if(is.na(match(file.ls[i],file.ls2))==F){#same file name exists in OutSpecPriors so alternative estimate available
			d2=read.table(file=paste0(OutDir2,file.ls[i]),header=T)
			irow=which(rownames(d2)=="Ntot")
			mu2[i]=d2$mean[irow]*yscale
			cl2[i,]=c(d2$X2.5.[irow],d2$X97.5.[irow])*yscale
		
		}
		StPos=as.integer(regexpr(pattern="_sum",file.ls[i]))-4
		yr[i]=substr(x=file.ls[i],start=StPos,stop=StPos+3)
		print(c(i,yr[i],d$Rhat[irow]))
	}
	
	if(ShowMissYrs==F){
		x=barplot(height=mu, names.arg=yr ,ylim=c(0,max(cl[,2])),bty='n',cex.names=0.9, border=NA,las=3,axis.lty=1,ylab="",xlab="",main=doTrib[itrib])
		for(i in 1:Nyrs) arrows(x0=x[i],x1=x[i],y0=cl[i,1],y1=cl[i,2],angle=90,code=3,length=0.025,col="black")
		if(ShowSpec==T){
			barplot(height=mu2, col=rgb(1,0,0,0.5),add=T,axes=F,border=NA)
			for(i in 1:Nyrs) arrows(x0=x[i],x1=x[i],y0=cl2[i,1],y1=cl2[i,2],angle=90,code=3,length=0.025,col="red")
		}
	
	} else {

		#This version includes blank years that were not estimated
		AllYr=seq(min(yr),max(yr))
		Nallyrs=length(AllYr)
		amu=vector(length=Nallyrs);amu=0;acl=matrix(nrow=Nallyrs,ncol=2,data=0)
		amu2=amu;acl2=acl
		for(i in 1:Nallyrs){
			ii=length(which(yr==AllYr[i]))
			if(ii>0){
				j=which(yr==AllYr[i])
				amu[i]=mu[j];acl[i,]=cl[j,]
				if(is.na(mu2[j])==F){
					amu2[i]=mu2[j]
					acl2[i,]=cl2[j,]
				}
			}
		}
		x=barplot(height=amu, names.arg=AllYr ,ylim=c(0,max(acl[,2])),bty='n',cex.names=0.9, las=3,axis.lty=1,ylab="",xlab="",main=doTrib[itrib])
		for(i in 1:Nallyrs) arrows(x0=x[i],x1=x[i],y0=acl[i,1],y1=acl[i,2],angle=90,code=3,length=0.025,col="black")
		if(ShowSpec==T){
			barplot(height=amu2, col=rgb(1,0,0,0.5),add=T,axes=F,border=NA)
			for(i in 1:Nallyrs) arrows(x0=x[i],x1=x[i],y0=acl2[i,1],y1=acl2[i,2],angle=90,code=3,length=0.025,col="red")
		}
		
	}
	abline(h=mean(mu,na.rm=T),lty=2,col="black")
	
	if(itrib==4 | itrib==8){
		mtext("Run Year (Nov-Dec yr 't-1' + Jan-May yr 't')",side = 1,line = -1, outer = T,cex=1.2,font=2)
		mtext("Abundance ('000s)",side=2, las=3,outer=T,,cex=1.2,font=2,line=-1)
	}

}