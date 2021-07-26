
library(FME)
mu <- 10
std <- 1
Nfun <- function(p) -2*log(dnorm(p, mean = mu, sd = std))


MCMC <- modMCMC (f = Nfun, p = 9.5, niter = 2000, jump = 5)
MCMC2 <- modMCMC (f = Nfun, p = 9.5, lower = 9, niter = 2000, jump = 5, updatecov = 10)
pri <- function(p) -2*log(dnorm(p, 8, 1))
MCMC3 <- modMCMC (f = Nfun, p = 9.5, niter = 2000, jump = 5, updatecov = 10, prior = pri)

par(mfrow = c(2,2))
hist(MCMC$pars, xlab="x", freq = FALSE, main = "unconstrained", xlim = c(6, 14))
hist(MCMC2$pars, xlab="x", freq = FALSE, main = "x>9", xlim = c(6, 14))
hist(MCMC3$pars, xlab="x", freq = FALSE, main = "pri(x)~N(8,1)", xlim = c(6, 14))
plot(MCMC3, mfrow = NULL, main = "AM")
mtext(outer = TRUE, line = -1.5, "N(10,1)", cex = 1.25)