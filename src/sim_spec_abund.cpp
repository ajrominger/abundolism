#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <random>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
Rcpp::List sim_spec_abund(const arma::vec& la,
                          const arma::vec& mu,
                          const arma::vec& g,
                          const arma::vec& m_prop,
                          const arma::vec& nu,
                          const arma::vec& tau,
                          const arma::vec& xi,
                          int np,
                          int nstep,
                          bool scale_setback) {

    // calculate dispersal rates
    vec m = m_prop % g;

    // number of simulations
    int nsim = la.n_elem;

    // output list
    Rcpp::List out(nsim);

    // total number of event types:
    // emigration (np long) + birth (np) + death (np) + speciation (np) +
    // immigration (1 long)
    int total_events = 4 * np + 1;

    // pre-allocate objects that will be updated in simulations
    uvec x(np);     // vector of local population sizes
    vec s(np);      // vector tracking time since incipient speciation
    vec stau(np);   // vector tracking time to wait for full speciation
    uvec sfull(np); // vector of booleans to indicate full speciation
    double tt;      // track total simulation time
    bool any_full_spec = false; // speciation tracker

    // set-up random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> sample_pop(0, np - 1);

    // loop over replicate simulations
    for(int r = 0; r < nsim; r++) {
        // re-set book keeping objects
        x.zeros();
        s.zeros();
        stau.zeros();
        sfull.zeros();
        tt = 0.0;
        any_full_spec = false;

        // pre-allocate this sim's matrix at max dim possible (nstep rows),
        // will trim to actual size before storing in the list
        mat sim_output(nstep, 4);

        // track this separate from `i` for trimming outside loop, relevant
        // in cases of early stopping
        int step_count = 0;

        for(int i = 0; i < nstep; i++) {
            // calculate rates
            vec this_la = conv_to<vec>::from(x) * la(r);
            vec this_mu = conv_to<vec>::from(x) * mu(r);
            vec this_m = conv_to<vec>::from(x) * m(r);
            vec this_nu = conv_to<vec>::from(x) * nu(r);
            double this_g = g(r);

            // create vector of probabilities of each event for
            // eventual sampling
            std::vector<double> probs(total_events);

            // fill probabilities: emigration, birth, death,
            //                     speciation, immigration
            for(int j = 0; j < np; j++) {
                probs[j] = this_m(j);         // emigration from pop j
                probs[np + j] = this_la(j);   // birth in pop j
                probs[2*np + j] = this_mu(j); // death in pop j
                probs[3*np + j] = this_nu(j); // speciation in pop j
            }

            probs[4*np] = this_g;             // immigration

            // sample event type using weighted discrete distribution
            std::discrete_distribution<int> sample_event(probs.begin(),
                                                         probs.end());

            int event_idx = sample_event(gen);

            // determine event type and population
            int e_type;     // will hold int representing event type
            int e_pop = -1; // -1 means immigration (updated below as need)

            // note on event coding: the first `np` elements of `event_idx`
            // are emigration, the next `np` elements are birth, etc
            // within each chunk of `np` elements, the relative index is
            // equal to local pop ID
            if(event_idx < np) {
                // emigration
                e_type = 0;
                e_pop = event_idx;
            } else if(event_idx < 2 * np) {
                // birth
                e_type = 1;
                e_pop = event_idx - np;
            } else if(event_idx < 3 * np) {
                // death
                e_type = 2;
                e_pop = event_idx - 2 * np;
            } else if(event_idx < 4 * np) {
                // speciation
                e_type = 3;
                e_pop = event_idx - 3 * np;
            } else {
                // immigration
                e_type = 4;
            }

            // sample sojourn time
            double total_rate = std::accumulate(probs.begin(),
                                                probs.end(),
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
            if(e_type == 0 || e_type == 4) { // emigration or immigration
                // choose receiving pop uniform random
                int receiving_pop = sample_pop(gen);

                // increment pop size
                x(receiving_pop) ++;

                // update wait time to speciation
                if(s(receiving_pop) > 0) {
                    if(scale_setback) {
                        stau(receiving_pop) += xi(r) / x(receiving_pop);
                    } else {
                        stau(receiving_pop) += xi(r);
                    }

                }
            } else if(e_type == 1) { // birth
                x(e_pop) ++;
            } else if(e_type == 2) { // death
                x(e_pop) --;
            } else if(e_type == 3) { // incipient speciation
                // only update time since incipient speciation if hasn't
                // already happened
                if(s(e_pop) == 0) {
                    s(e_pop) = st;
                }
            }

            // update total simulation time
            tt += st;

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

            // record this step
            sim_output(step_count, 0) = r;
            sim_output(step_count, 1) = sum(conv_to<vec>::from(x));
            sim_output(step_count, 2) = st;
            sim_output(step_count, 3) = any_full_spec ? 1.0 : 0.0;
            step_count++;


            if(any_full_spec) {
                break;
            }
        }

        // trim to actual steps taken and store
        out[r] = sim_output.rows(0, step_count - 1);
    }

    return out;
}
