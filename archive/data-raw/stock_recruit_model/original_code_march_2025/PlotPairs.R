
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
  r <- (cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = 1.3)#cex.cor * r
}


par(mfcol=c(1,1),mai=c(0.85,0.85,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))

dp=as.data.frame(fit, pars = c("alpha","beta","gamma","sd_pro"))

icol=which(names(dp)=="alpha");alpha=dp[,icol]
icol=which(names(dp)=="beta");beta=dp[,icol]
icol=which(names(dp)=="gamma");gamma=dp[,icol]
icol=which(names(dp)=="sd_pro");sd_pro=dp[,icol]
xpars=cbind(alpha,beta,gamma,sd_pro)

xlabels=c(expression(alpha),expression(beta),expression(gamma),expression(paste(sigma[p])))
pairs(x=xpars,cex=0.5,labels=xlabels,diag.panel=panel.hist,upper.panel=panel.cor)