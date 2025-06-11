#' @title r_sim
#' @description
#' A pure R implementation to simulate speciation in a meta-population
#' birth-death-immigration model; mostly for testing against rcpp version
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
#' @returns a matrix

r_sim <- function(la, mu, g, m_prop, nu, tau, xi, np, nstep) {
    # dispersal rate between communities
    m <- m_prop * g

    # vector of events for sampling during simulation
    e <- c(paste("emm", 1:np, sep = "_"),
           paste("b", 1:np, sep = "_"),
           paste("d", 1:np, sep = "_"),
           paste("spec", 1:np, sep = "_"),
           "imm")

    # matrix to store output
    nsim <- length(la)
    output <- matrix(NA, nrow = nsim, ncol = 3)

    # loop over replicate simulations
    for(r in 1:nsim) {
        # vector of local population sizes
        x <- rep(0, np)

        # running sum (for calculating mean abund) across all pops
        xmean <- sum(x)

        # vector traking time since incipient speciation
        s <- numeric(np)

        # vector traking time gotta wait to full speciation
        stau <- rep(tau[r], np)

        # vector of booleans if there has been full speciation
        sfull <- rep(FALSE, np)

        # keep track of total sim time
        tt <- 0

        # loop over steps to simulate process
        for(i in 2:nstep) {
            # rates
            this_la <- x * la[r]
            this_mu <- x * mu[r]
            this_m <- x * m[r]
            this_nu <- x * nu[r]
            this_g <- g[r]

            # sample event
            this_e <- sample(e, 1, prob = c(this_m, this_la, this_mu,
                                            this_nu, this_g))

            # parse event
            e_type <- gsub("_.*", "", this_e) # what kind of event

            if(e_type != "imm") {
                e_pop <- as.numeric(gsub(".*_", "", this_e)) # who originates the event
            }

            # sample sojourn time
            st <- rexp(1, this_m + this_la + this_mu + this_nu + g)

            # if there has been incipient speciation, add st to time since
            # incipient speciation
            s[s > 0] <- s[s > 0] + st

            # update based on event type
            if(e_type == "emm" | e_type == "imm") {
                receiving_pop <- sample(np, 1)

                # add to pop size
                x[receiving_pop] <- x[receiving_pop] + 1

                # update wait time to speciation
                if(s[receiving_pop] > 0) {
                    stau[receiving_pop] <- stau[receiving_pop] +
                        xi[r] / x[receiving_pop]
                }

            } else if(e_type == "b") {
                x[e_pop] <- x[e_pop] + 1

            } else if(e_type == "d") {
                x[e_pop] <- x[e_pop] - 1

            } else { # e_type == "spec"
                # only update time since incipient speciation if hasnʻt alrady happened
                if(s[e_pop] == 0) {
                    s[e_pop] <- st
                }
            }

            # update total sim time
            tt <- tt + st

            # update running population sum
            xmean <- xmean + sum(x)

            # if local extirpation, reset speciation clock
            stau[x == 0] <- tau[r]
            s[x == 0] <- 0

            # check for full speciation
            if(any(s > stau)) {
                sfull[s > stau] <- TRUE

                break
            }
        }

        # record this simulation
        output[r, ] <- c(tt, xmean / i / np, 1 * any(sfull))
    }

    colnames(output) <- c("time", "xmean", "full_spec")

    return(output)
}
