#' @title sim_BDI_spec
#' @description
#' A fast simulation of speciation in a meta-population
#' birth-death-immigration model
#' @param la vector of bird rates
#' @param mu vector of death rates
#' @param g vector of immigration rates from global population
#' @param m_prop vector of proportions (local immigration):(global immigartion)
#' @param nu vector of rates of incipient speciation
#' @param tau waiting time to full speciation
#' @param xi vector of weights given to new immigrants in slowing progress to
#'           full speciation
#' @param np number of local populations in the meta population
#' @param nstep scalar, number of simulation steps
#'
#' @returns a data.frame with columns for input parameter values and output
#' results
#' - `la`
#' - `mu`
#' - `g`
#' - `m_prop`
#' - `nu`
#' - `tau`
#' - `xi`
#' - `np`
#' - `nstep`
#' - `time`: simulation time step, if speciation occurred this number will
#' be `< nstep`, if not this number will be `== nstep`
#' - `mean_pop_size`: the mean population size across all time steps and
#' populations
#' - `speciation`: 0 (no speciation) or 1 (yes specaition)
#' @md
#'
#' @export

sim_BDI_spec <- function(la, mu, g, m_prop, nu, tau, xi, np, nstep) {
    sim <- sim_spec_abund(la, mu, g, m_prop, nu, tau, xi, np, nstep)

    sim <- as.data.frame(sim)

    names(sim) <- c("time", "niter", "mean_pop_size", "speciation")

    pars <- data.frame(la, mu, g, m_prop, nu, tau, xi, np, nstep)

    sim <- cbind(pars, sim)

    return(sim)
}

