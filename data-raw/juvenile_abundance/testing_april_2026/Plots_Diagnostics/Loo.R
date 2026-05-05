library("loo") #https://mc-stan.org/loo/articles/loo2-with-rstan.html


#Do LOO calculations
log_lik <- extract_log_lik(pcap, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
loo <- loo(log_lik, r_eff = r_eff)
#looic_stats[irun,]=loo$estimates[3,]#save the looic estimate and its SE
#loo$pointwise #yr-specific looic values highlighting influential points
print(loo)
