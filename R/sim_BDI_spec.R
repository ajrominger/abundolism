#' @title Simulate a birth-death-speciation model
#' @description
#' A fast simulation of speciation in a meta-population
#' birth-death-immigration model
#' @param la vector of bird rates
#' @param mu vector of death rates
#' @param g vector of immigration rates from global population
#' @param m_prop vector of proportions (local immigration):(global immigration)
#' @param nu vector of rates of incipient speciation
#' @param tau waiting time to full speciation
#' @param xi vector of weights given to new immigrants in slowing progress to
#'           full speciation
#' @param np number of local populations in the meta population
#' @param nstep scalar, number of simulation steps
#' @param w scalar, number of simulation steps in window for sliding average
#'
#' @returns a data.frame with columns for input parameter values and output
#' results
#' - `sim_id`
#' - `la`
#' - `mu`
#' - `g`
#' - `m_prop`
#' - `nu`
#' - `tau`
#' - `xi`
#' - `np`
#' - `nstep`
#' - `w`
# - `time`: simulation time step, if speciation occurred this number will
#           be `< nstep`, if not this number will be `== nstep`
#' - `mean_pop_size`: the mean population size across all time steps and
#'                    populations
#' - `speciation`: 0 (no speciation) or 1 (yes speciation)
#'
#' @md
#'
#' @export

sim_BDI_spec <- function(la, mu, g, m_prop, nu, tau, xi, np, nstep, w) {
    # if(length(nstep) == 1) nstep <- rep(nstep, length(la))

    raw_sim <- sim_spec_abund(la, mu, g, m_prop, nu, tau, xi, np, nstep)

    sim <- lapply(1:length(raw_sim), function(i) {
        # this sim
        x <- raw_sim[[i]]

        # breakpoints across windows
        b <- seq(from = nrow(x), to = 0, by = -w)
        l <- split(x[, 2], cut(1:nrow(x), b))

        # pop averages and yes/no speciation in windows
        pop_m <- vapply(l, mean, FUN.VALUE = 0)
        spec <- rep(0, length(pop_m))
        if(any(x[, 4] > 0)) spec[length(pop_m)] <- 1

        data.frame(
            sim_id = x[1, 1] + 1, # +1 for r-style index
            la = la[i], mu = mu[i], g = g[i],
            m_prop = m_prop[i], nu = nu[i],
            tau = tau[i], xi = xi[i],
            np = np[i], nstep = nstep, w = w,
            mean_pop_size = pop_m,
            speciation = spec
        )
    }) |>
        do.call(rbind, args = _)

    rownames(sim) <- NULL

    return(sim)
}

