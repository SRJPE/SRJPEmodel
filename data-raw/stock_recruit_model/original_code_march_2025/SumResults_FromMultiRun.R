#Calculate fit and loo ic and key parameter summaries from model that was just run

#function to calculated overlap in two LOO IC normal distributions
min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
  f1 <- dnorm(x, mean=mu1, sd=sd1)
  f2 <- dnorm(x, mean=mu2, sd=sd2)
  pmin(f1, f2)
}


if(CoVarNm!="mu" & CoVarNm!="null"){
  dp=as.data.frame(fit, pars = c("gamma")) #for discrete models this only shows gamma[1]
  maxgcol=dim(dp)[2]
  gamma_stats[irun,1:3]=round(as.double(quantile(dp[,maxgcol],probs=c(0.025,0.5,0.975))),digits=2)
  gamma_stats[irun,4]=round(abs(sd(dp[,maxgcol])/mean(dp[,maxgcol])),digits=2)
} else {
  gamma_stats[irun,1:4]=rep(NA,4)
}


obslgRS=log(R/SP);  n[irun]=length(obslgRS)
if(CoVarNm!="mu" & CoVarNm!="null"){
  dp=as.data.frame(fit, pars = c("alpha","beta","gamma"))
  predlgRS=matrix(nrow=dim(dp)[1],ncol=Nyrs)
  
  for(iyr in 1:Nyrs){
    
    if(CoVarNm!="WY2" & CoVarNm!="WY3") {
      predlgRS[,iyr]=dp$alpha + dp$beta*SP[iyr] + dp$gamma*X[iyr]
      
    } else { #wy2 or wy3
      
      if(X[iyr]==0){
        predlgRS[,iyr]=dp$alpha + dp$beta*SP[iyr]
      } else {
        if(Ndisc>1){
          vnm=paste0("gamma[",X[iyr],"]")
        } else {
          vnm="gamma"
        }
        icol=which(names(dp)==vnm)
        predlgRS[,iyr]=dp$alpha + dp$beta*SP[iyr] + dp[,icol]
      }
    }
  }
  r2[irun]=round(cor(obslgRS,colMeans(predlgRS))^2,digits=2)

} else if (CoVarNm=="null"){
    dp=as.data.frame(fit, pars = c("alpha","beta"))
    predlgRS=matrix(nrow=dim(dp)[1],ncol=Nyrs)
    for(iyr in 1:Nyrs){
      predlgRS[,iyr]=dp$alpha + dp$beta*SP[iyr]
    }
    r2[irun]=round(cor(obslgRS,colMeans(predlgRS))^2,digits=2)

} else { #mean model
  r2[irun]=0
  #dp=as.data.frame(fit, pars = c("alpha"))
  #predlgRS=matrix(nrow=dim(dp)[1],ncol=Nyrs)
  #for(iyr in 1:Nyrs){
  #  predlgRS[,iyr]=dp$alpha*rnorm(n=dim(dp)[1],mean=1,sd=0.001) #fuzzy it up so cor doesn't bomb (happens if no variation in predlgRS)
  #}
  #r2[irun]=round(cor(obslgRS,colMeans(predlgRS))^2,digits=2)
}

dp=as.data.frame(fit,pars=c("sd_pro"))
sd_pro[irun]=mean(dp[,1])

#Do LOO calculations
log_lik <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
loo <- loo(log_lik, r_eff = r_eff)
looic_stats[irun,]=loo$estimates[3,]#save the looic estimate and its SE
#loo$pointwise #yr-specific looic values highlighting influential points

if(irun==Ncovars){
  im=which.min(looic_stats[,1]) #index for lowest loo ic model
  mu1=looic_stats[im,1];sd1=looic_stats[im,2]
  for(jrun in 1:Ncovars){
    looic_stats[jrun,3]=looic_stats[jrun,1]-looic_stats[im,1]
    
    #pover_stats=integrate(min.f1f2, -Inf, Inf, mu1=looic_stats[im,1], mu2=looic_stats[jrun,1], sd1=looic_stats[im,2], sd2=looic_stats[jrun,2])
    #Need narrower range for integration or it evalutes to zero when se's of LOO IC are small (knights landing)
    pover_stats=integrate(min.f1f2, -2*looic_stats[im,1], 2*looic_stats[im,1], mu1=looic_stats[im,1], mu2=looic_stats[jrun,1], sd1=looic_stats[im,2], sd2=looic_stats[jrun,2])
    looic_stats[jrun,4]=round(as.double(pover_stats[1]),digits=3)
  }
  
  fn=paste0(OutDir,DoSite,"_",Am,"_","Sum.dat")
  write.table(file=fn,cbind(unCoVarNm,n,r2,sd_pro,gamma_stats,looic_stats),row.names=F,col.names=c("Covariate","Nyrs","r2","sd_pro","LCL","Median","UCL","CV","mu_looic","se_looic","delta_loo","pover_loo"))
}
