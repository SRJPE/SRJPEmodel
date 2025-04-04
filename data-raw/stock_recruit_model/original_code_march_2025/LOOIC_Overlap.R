min.f1f2 <- function(x, mu1, mu2, sd1, sd2) {
    f1 <- dnorm(x, mean=mu1, sd=sd1)
    f2 <- dnorm(x, mean=mu2, sd=sd2)
    pmin(f1, f2)
}

#mu1 <- 2;sd1 <- 2;mu2 <- 1;sd2 <- 1
mu1=21.46;sd1=4.02

#mu2=22.52;sd2=6.01
mu2=26.63;sd2=4.64


xs <- seq(min(mu1 - 3*sd1, mu2 - 3*sd2), max(mu1 + 3*sd1, mu2 + 3*sd2), .01)
f1 <- dnorm(xs, mean=mu1, sd=sd1)
f2 <- dnorm(xs, mean=mu2, sd=sd2)

plot(xs, f1, type="l", ylim=c(0, max(f1,f2)), ylab="density")
lines(xs, f2, lty="dotted")
ys <- min.f1f2(xs, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)
xs <- c(xs, xs[1])
ys <- c(ys, ys[1])
polygon(xs, ys, col="gray")

### only works for sd1 = sd2
SMD <- (mu1-mu2)/sd1
2 * pnorm(-abs(SMD)/2)

### this works in general
poverlap=integrate(min.f1f2, -Inf, Inf, mu1=mu1, mu2=mu2, sd1=sd1, sd2=sd2)
print(poverlap)
