model{

  # Prior distributions of the model parameters
  tau.p ~dgamma(0.001,0.001)
  mu.p~dbeta(1,1)I(0.001,)
  L.mu.p <- logit(mu.p)

  mu.kappa ~ dgamma(1,0.001)
  tau.kappa ~ dgamma(0.001,0.001)I(0.001,)
  alpha <- mu.kappa * beta
  beta <- mu.kappa * tau.kappa

  for (i in 1:3){
    mu.A[i]~dunif(0,500000)
    L.mu.A[i]<-log(mu.A[i])
    tau.A[i]~dgamma(0.001,0.001)
  }

  # calibration between adult counts and redd counts
  # T = 8 years, 2001 to 2008; i = rivers (1=Faughan; 2=Finn; 3=Roe)
  for (t in 1:T){
    for (i in 1:3){
      L.p[t,i] ~ dnorm(L.mu.p, tau.p)
      logit(p[t,i]) <- L.p[t,i]
      A[t,i] ~ dlnorm(L.mu.A[i],tau.A[i])
      S[t,i] <- A[t,i] - C[t,i]
      kappa[t,i] ~ dgamma(alpha, beta)I(0.001,)
      lambda[t,i] <- (A[t,i] - C[t,i]) * p[t,i] * kappa[t,i]
      R[t,i] ~dpois(lambda[t,i])
    }
  }
  # time series
  # Y = 42 years, 1959 to 2000; i = rivers (1=Faughan; 2=Finn; 3=Roe)
  for (y in 1:Y){
    for (i in 1:3){
      L.p.old[y,i] ~ dnorm(L.mu.p, tau.p)
      logit(p.old[y,i]) <- L.p.old[y,i]

      A.old[y,i] ~ dlnorm(L.mu.A[i],tau.A[i])

      kappa.old[y,i] ~ dgamma(alpha, beta)I(0.001,)
      lambda.redds.old[y,i] <- (A.old[y,i] – C[t,i]) * p.old[y,i] * kappa.old[y,i]
      Redds.old[y,i] ~dpois(lambda.redds.old[y,i])
    }

  }
}
