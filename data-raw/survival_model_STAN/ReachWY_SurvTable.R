#Calclualte S100 by reach and water year type

library("rstan")

inv_logit <-function(x){
	return(exp(x)/(1+exp(x)))
}


varnm1="S_bCov";varnm2="RE_sd"

for(irun in 1:7){
  ModNm=switch(irun,"NoCov","CovWY2","CovWY2_Reach","CovWY3_Reach","CovInd_Reach","CovInd_Reach","CovInd_Reach")

  fnext=""
  if(irun==5){
    fnext="_Q"
  } else if (irun==6){
    fnext="_Vel"
  } else if (irun==7){
    fnext="_Wt"
  }
  
  fnNm=paste0("Results/fit_",ModNm,fnext,".Rdata")
  load(file=fnNm)
  
  S_bReach=as.data.frame(fit,pars=c("S_bReach"))
  Nsim=dim(S_bReach)[1]
  S_bR=matrix(nrow=Nsim,ncol=3)
  for(j in 1:3){
      vname=paste0("S_bReach[",j,"]")
      icol=which(names(S_bReach)==vname)
      S_bR[,j]=S_bReach[,icol]
  }
  
  if(irun==1){
   	muS100=vector(length=3);sdS100=muS100
   	for(j in 1:3){
   		muS100[j]=round(mean(inv_logit(S_bR[,j])),digits=2)
   		sdS100[j]=round(sd(inv_logit(S_bR[,j])),digits=2)
   	}
 	
  } else {
  
  	S_bCov=as.data.frame(fit, pars = c(varnm1))  
  	
  	if(irun==2){
  		S_bC=unlist(S_bCov)
  		muS100=matrix(nrow=2,ncol=3);sdS100=muS100
  		for(j in 1:3){
  			muS100[1,j]=round(mean(inv_logit(S_bR[,j])),digits=2) #dry
  			sdS100[1,j]=round(sd(inv_logit(S_bR[,j])),digits=2)
  			muS100[2,j]=round(mean(inv_logit(S_bR[,j]+S_bC[])),digits=2) #wet
  			sdS100[2,j]=round(sd(inv_logit(S_bR[,j]+S_bC[])),digits=2)
  		}
    	
  	} else if (irun==3){
	    muS100=matrix(nrow=2,ncol=3);sdS100=muS100
	    for(j in 1:3){
	      vname=paste0(varnm1,"[",j,"]")
	      icol=which(names(S_bCov)==vname)
	      S_bC=S_bCov[,icol]
	      muS100[1,j]=round(mean(inv_logit(S_bR[,j])),digits=2) #dry
	      sdS100[1,j]=round(sd(inv_logit(S_bR[,j])),digits=2)
	      muS100[2,j]=round(mean(inv_logit(S_bR[,j]+S_bC[])),digits=2) #wet
	      sdS100[2,j]=round(sd(inv_logit(S_bR[,j]+S_bC[])),digits=2)
	    }
	
	} else if (irun==4){
	    muS100=matrix(nrow=3,ncol=3);sdS100=muS100
	    for(j in 1:3){
	    	muS100[1,j]=round(mean(inv_logit(S_bR[,j])),digits=2) #dry
		    sdS100[1,j]=round(sd(inv_logit(S_bR[,j])),digits=2)
	    	
	    	for(iwy in 1:2){
	         	vname=paste0(varnm1,"[",iwy,",",j,"]")
	        	icol=which(names(S_bCov)==vname)
	        	S_bC=S_bCov[,icol]
	        	muS100[iwy+1,j]=round(mean(inv_logit(S_bR[,j]+S_bC[])),digits=2) #wet
	      		sdS100[iwy+1,j]=round(sd(inv_logit(S_bR[,j]+S_bC[])),digits=2)
	       }#next ity
	    }#next j
  
  	} else if (irun>=5){ #show survival at 1SD below, 0, and 1SD above mean covariate value
	    muS100=matrix(nrow=3,ncol=3);sdS100=muS100
	    S_bC=unlist(S_bCov)
	    
	    for(ix in 1:3){
	    	X=switch(ix,-1,0,1)
	    	for(j in 1:3){
	       	  muS100[ix,j]=round(mean(inv_logit(S_bR[,j]+S_bC[]*X)),digits=2) #dry
	    		  sdS100[ix,j]=round(sd(inv_logit(S_bR[,j]+S_bC[]*X)),digits=2)
	       }#next j
	    }#next i
	    
  	}#irun if statement
   	
   } #if irun>1
 
  RE_sd=unlist(as.data.frame(fit, pars = c(varnm2)))
  mu_RE_sd=round(mean(RE_sd),digits=2); sd_RE_sd=round(sd(RE_sd),digits=2)
  
  print("")
  print(ModNm)
  print(muS100);print(sdS100)
  print(c(mu_RE_sd,sd_RE_sd))
  
} #next irun
