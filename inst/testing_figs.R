library(abundolism)
library(ggplot2)
library(ggpointdensity)

nrep <- 500
sim <- sim_BDI_spec(la = runif(nrep, 0.1, 4), mu = runif(nrep, 0.1, 4),
                    g = runif(nrep, 0.1, 1), m_prop = runif(nrep, 0, 0.1),
                    nu = runif(nrep, 0, 0.01), #nu = runif(nrep, 0.01, 0.01),
                    tau = runif(nrep, 20, 200), xi = runif(nrep, 10, 100),
                    np = 2, nstep = 20000) #20000

ggplot(sim, aes(mu / la, mean_pop_size)) +
    geom_pointdensity(method = "kde2d") +
    scale_color_viridis_c() +
    facet_wrap(vars(speciation)) +
    scale_x_log10() +
    scale_y_log10()

ggplot(sim, aes(time, mean_pop_size)) +
    geom_pointdensity(method = "kde2d") +
    scale_color_viridis_c() +
    facet_wrap(vars(speciation)) +
    scale_x_log10() +
    scale_y_log10()

ggplot(sim, aes(mu / la)) +
    geom_histogram() +
    scale_x_log10() +
    facet_grid(rows = vars(speciation))


ggplot(sim, aes(mean_pop_size, speciation)) +
    geom_pointdensity(method = "kde2d") +
    scale_color_viridis_c() +
    scale_x_log10()


sim <- sim_BDI_spec(la = rep(1, nrep), mu = rep(1, nrep),
                    g = rep(1, nrep), m_prop = rep(0.1, nrep),
                    nu = rep(0.001, nrep), tau = rep(20, nrep), xi = rep(4, nrep),
                    np = 2, nstep = 20000)
