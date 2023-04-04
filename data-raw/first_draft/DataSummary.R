#Summarize the RST data

rm(list=ls(all=TRUE)) 

fn="RST_Input.txt"	#all stream_sites available as of Aug 02

Fwk=45; Lwk=22
Week=c(seq(1,Lwk),seq(Fwk,53));Nwks=length(Week)

SampProp=0.75

d=read.table(file=fn,header=T)
Tribs=unique(d$Trib)
Ntribs=length(Tribs)

Nyr=vector(length=Ntribs);Fyr=Nyr;Lyr=Nyr;PassYrs=Nyr;Nmr=Nyr
Nyrmr=matrix(data=NA,nrow=30,ncol=Ntribs)
MRstat=matrix(nrow=Ntribs,ncol=3)

for(i in 1:Ntribs){
	d1=subset(d,Trib==Tribs[i])
	
	d1a=subset(d1,is.na(Rel1)==F)
	Nmr[i]=dim(d1a)[1]
	
	unq_yr=unique(d1$RunYr)
	Nyr[i]=length(unq_yr);Fyr[i]=min(unq_yr);Lyr[i]=max(unq_yr)
	
	k=0
	for(j in 1:Nyr[i]){
		d2=subset(d1,RunYr==unq_yr[j] & is.na(ua)==F)
		if(dim(d2)[1]/Nwks>=SampProp) k=k+1
	
		d3=subset(d1,RunYr==unq_yr[j] & is.na(Rel1)==F)
		Nyrmr[j,i]=dim(d3)[1]
	}
	PassYrs[i]=k	#round(k/Nyr[i],digits=2)
	
	MRstat[i,]=round(as.double(quantile(Nyrmr[,i],probs=c(0.2,0.5,0.8),na.rm=T)),digits=1)
}
out=data.frame(cbind(Tribs,Fyr,Lyr,Nyr,PassYrs,Nmr,MRstat))
names(out)=c("Str_Site","Fyr","Lyr","Nyrs","PassYrs","Nmr","MRwksq02","MRwksq05","MRwksq08")
write.table("DataSum.out",x=out)
