library("rstan")

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T);par(las=1)


par(mfcol=c(1,1),mai=c(0.85,0.85,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))

panel.hist <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.3)#cex.cor * r
}

dp=as.data.frame(fit, pars = c("phi","lambda","sd_pro","bCov","muRT","sigmaRT","rho"))#,"rho_pro"))

icols=c(which(names(dp)=="muRT[1]"),icols=which(names(dp)=="muRT[2]"));muRT=dp[,icols]
icols=c(which(names(dp)=="sigmaRT[1]"),icols=which(names(dp)=="sigmaRT[2]"));sigmaRT=dp[,icols]
icol=which(names(dp)=="rho");rho=dp[,icol]
icol=which(names(dp)=="sd_pro");sd_pro=dp[,icol]
#icol=which(names(dp)=="rho_pro");rho_pro=dp[,icol]
icol=which(names(dp)=="bCov[1]");bCov=dp[,icol]
for(iyr in 1:1){#1:Nyrs
  icol=which(names(dp)==paste0("phi[",iyr,"]"));phi=dp[,icol]
  icol=which(names(dp)==paste0("lambda[",iyr,"]"));lambda=dp[,icol]
  #icol=which(names(dp)==paste0("sd_pro[",iyr,"]"));sd_pro=dp[,icol]
  
  #xpars=cbind(phi,lambda,sd_pro,bCov, rho)#, rho_pro)#muRT,sigmaRT,rho)
  xpars=cbind(muRT[1],sigmaRT,rho,bCov)#, rho_pro)#muRT,sigmaRT,rho)
  
  pairs(x=xpars,pch=19,cex=0.2,diag.panel=panel.hist,upper.panel=panel.cor,main=year[iyr])
}
