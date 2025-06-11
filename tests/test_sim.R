library(abundolism)

la = c(1, 10)
mu = c(3, 3)
g = c(2, 3)
m_prop = c(0.1, 0.1)
nu = c(0.1, 0.1)
tau = c(1, 10)
xi = c(1, 1)
np = 2
nstep = 400

nt <- proc.time()
x <- abundolism:::sim_spec_abund(la = la, mu = mu, g = g, m_prop = m_prop,
                                 nu = nu, tau = tau, xi = xi, np = np,
                                 nstep = nstep)
proc.time() - nt
x

nt <- proc.time()
y <- abundolism:::r_sim(la = la, mu = mu, g = g, m_prop = m_prop,
                        nu = nu, tau = tau, xi = xi, np = np,
                        nstep = nstep)
proc.time() - nt

x
y




nrep <- 1000
la = runif(nrep, 0.1, 10)
mu = runif(nrep, 0.1, 10)
g = la * runif(nrep, 0, 0.1)
m_prop = runif(nrep, 0, 0.1)
nu = runif(nrep, 0, 0.1)
tau = 10 / (la + mu + runif(nrep, 0, 1))
xi = rep(1, nrep)
np = 2
nstep = 10000


x <- abundolism:::sim_spec_abund(la = la, mu = mu, g = g, m_prop = m_prop,
                                 nu = nu, tau = tau, xi = xi, np = np,
                                 nstep = nstep)



z <- sim_BDI_spec(la = la, mu = mu, g = g, m_prop = m_prop,
                  nu = nu, tau = tau, xi = xi, np = np,
                  nstep = nstep)

plot(z$mean_pop_size, z$speciation, log = "x")



