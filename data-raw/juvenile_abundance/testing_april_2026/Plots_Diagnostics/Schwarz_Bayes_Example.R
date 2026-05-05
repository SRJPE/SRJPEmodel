N.Years <- 10
N.Sites <- 5
N.Weeks <- 10
year.effect.sd <- 0
site.effect.sd <- 2
site.year.effect.sd <- 0
site.year.week.effect.sd <- 0


eff.data <- expand.grid(Year=1:N.Years, Site=1:N.Sites, Week=1:N.Weeks)
eff.data$logit.P.base <- log(.025/(1-.025))
year.effects <- data.frame(Year=1:N.Years, year.effect=rnorm(N.Years, mean=0, sd=year.effect.sd))
site.year.effects <- expand.grid(Year=1:N.Years, Site=1:N.Sites)
site.year.effects$site.year.effect = rnorm(N.Years*N.Sites, mean=0, sd=site.year.effect.sd)
site.effects <- data.frame(Site=1:N.Sites, site.effect=rnorm(N.Sites, mean=0, sd=site.effect.sd))
eff.data <- merge(eff.data, year.effects)
eff.data <- merge(eff.data, site.year.effects)
eff.data <- merge(eff.data, site.effects)

eff.data$logit.P <- eff.data$logit.P.base + eff.data$year.effect + eff.data$site.effect + eff.data$site.year.effect + rnorm(nrow(eff.data), mean=0, sd=site.year.week.effect.sd)
eff.data$P <- 1/(1+exp(-eff.data$logit.P)) #Caps

eff.data$n <- rep(c(seq(20,100,40), seq(150,1000,50)), length.out=nrow(eff.data)) #releases
eff.data$x <- rbinom(nrow(eff.data), size=eff.data$n, prob=eff.data$P) #simulated recaptures
