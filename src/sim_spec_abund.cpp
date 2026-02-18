#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <random>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat sim_spec_abund(const arma::vec& la,
                         const arma::vec& mu,
                         const arma::vec& g,
                         const arma::vec& m_prop,
                         const arma::vec& nu,
                         const arma::vec& tau,
                         const arma::vec& xi,
                         int np,
                         int nstep) {

    // calculate dispersal rates
    vec m = m_prop % g;  // element-wise multiplication

    // number of simulations
    int nsim = la.n_elem;

    // output matrix
    mat output(nsim, 3);
    // output.fill(datum::nan);

    // total number of event types:
    // emigration (np) + birth (np) + death (np) + speciation (np) +
    // immigration (1)
    int total_events = 4 * np + 1;

    // pre-allocate objects that will be updated in simulations
    uvec x(np);
    double xmean;
    vec s(np);
    vec stau(np);
    uvec sfull(np);
    double tt;
    bool any_full_spec = false;

    // set-up random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> sample_pop(0, np - 1);

    // loop over replicate simulations
    for(int r = 0; r < nsim; r++) {

        // re-set book keeping objects
        x.zeros();     // vector of local population sizes
        xmean = 0.0;   // running sum for calculating mean abundance
        s.zeros();     // vector tracking time since incipient speciation
        stau.zeros();  // vector tracking time to wait for full speciation
        sfull.zeros(); // vector of booleans to indicate full speciation
        tt = 0.0;      // track total simulation time
        any_full_spec = false; // speciation tracker

        // loop over simulation steps
        int i; // we are going to use this in mean abund calc
        for(i = 0; i < nstep; i++) {

            // calculate rates
            vec this_la = conv_to<vec>::from(x) * la(r);
            vec this_mu = conv_to<vec>::from(x) * mu(r);
            vec this_m = conv_to<vec>::from(x) * m(r);
            vec this_nu = conv_to<vec>::from(x) * nu(r);
            double this_g = g(r);

            // create probability vector for event sampling
            std::vector<double> probs(total_events);

            // fill probabilities: emigration, birth, death, speciation, immigration
            for(int j = 0; j < np; j++) {
                probs[j] = this_m(j);         // emigration from pop j
                probs[np + j] = this_la(j);   // birth in pop j
                probs[2*np + j] = this_mu(j); // death in pop j
                probs[3*np + j] = this_nu(j); // speciation in pop j
            }

            probs[4*np] = this_g;             // immigration

            vec arma_probs(probs.data(), probs.size(), false);

            // sample event type using weighted discrete distribution
            std::discrete_distribution<int> sample_event(probs.begin(),
                                                         probs.end());

            int event_idx = sample_event(gen);

            // determine event type and population
            int e_type;
            int e_pop = -1; // -1 is filler value for immigration

            if(event_idx < np) {
                e_type = 0;  // emigration
                e_pop = event_idx;
            } else if(event_idx < 2*np) {
                e_type = 1;  // birth
                e_pop = event_idx - np;
            } else if(event_idx < 3*np) {
                e_type = 2;  // death
                e_pop = event_idx - 2*np;
            } else if(event_idx < 4*np) {
                e_type = 3;  // speciation
                e_pop = event_idx - 3*np;
            } else {
                e_type = 4;  // immigration
            }

            // sample sojourn time
            double total_rate = std::accumulate(probs.begin(), probs.end(),
                                                0.0);
            std::exponential_distribution<double> exp_dist(total_rate);
            double st = exp_dist(gen);

            // update time since incipient speciation
            for(int j = 0; j < np; j++) {
                if(s(j) > 0) {
                    s(j) += st;
                }
            }

            // update pops based on event type
            if(e_type == 0 || e_type == 4) {  // emigration or immigration
                // choose receiving population randomly
                int receiving_pop = sample_pop(gen);

                // add to population size
                x(receiving_pop) ++;

                // update wait time to speciation
                if(s(receiving_pop) > 0) {
                    stau(receiving_pop) += xi(r) / x(receiving_pop);
                }

            } else if(e_type == 1) {  // birth
                x(e_pop) ++;

            } else if(e_type == 2) {  // death
                x(e_pop) --;

            } else if(e_type == 3) {  // incipient speciation
                // only update time since incipient speciation if hasn't
                // already happened
                if(s(e_pop) == 0) {
                    s(e_pop) = st;
                }
            }

            // update total simulation time
            tt += st;

            // update running population sum
            xmean += sum(conv_to<vec>::from(x));

            // if local extirpation, reset speciation clock
            for(int j = 0; j < np; j++) {
                if(x(j) == 0) {
                    stau(j) = tau(r);
                    s(j) = 0;
                }
            }

            // check for full speciation
            sfull = s > stau;
            any_full_spec = any(sfull);

            if(any_full_spec) {
                break;
            }
        }

        // record this simulation
        output(r, 0) = tt;                        // time
        output(r, 1) = xmean / i / np;            // mean abundance
        output(r, 2) = any_full_spec ? 1.0 : 0.0; // full speciation = 1
    }

    return output;
}
