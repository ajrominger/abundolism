#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <random>

using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat check_spec(const arma::vec& la,
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
    // int nsim = la.n_elem;

    // output matrix (different dim from main function)
    mat output(nstep, 3);


    // total number of event types:
    // emigration (np long) + birth (np) + death (np) + speciation (np) +
    // immigration (1 long)
    int total_events = 4 * np + 1;

    // pre-allocate objects that will be updated in simulations
    uvec x(np);     // vector of local population sizes
    // double xmean;   // running sum for calculate mean abund (across pops)
    vec s(np);      // vector tracking time since incipient speciation
    vec stau(np);   // vector tracking time to wait for full speciation
    uvec sfull(np); // vector of booleans to indicate full speciation
    // double tt;      // track total simulation time
    // bool any_full_spec = false; // speciation tracker

    // set-up random number generators
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> sample_pop(0, np - 1);

    // in full function we loop over r
    int r = 0;

    // re-set book keeping objects
    // x.zeros();
    x.fill(10); // fill with abund > 0 so we have chance of all probs
    // xmean = 0.0;
    s.zeros();
    stau.zeros();
    sfull.zeros();
    // tt = 0.0;
    // any_full_spec = false;

    // loop over simulation steps
    int i; // we are going to use this in mean abund calc
    for(i = 0; i < nstep; i++) {

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

        // not used
        // vec arma_probs(probs.data(), probs.size(), false);

        // sample event type using weighted discrete distribution
        std::discrete_distribution<int> sample_event(probs.begin(),
                                                     probs.end());

        int event_idx = sample_event(gen);

        // determine event type and population
        int e_type;
        int e_pop = -1; // -1 means immigration (updated below as need)

        if(event_idx < np) { // is this right, or should be np - 1?
            // emigration
            e_type = 0;
            e_pop = event_idx; // `event_idx` set up so pop id = f(idx)
        } else if(event_idx < 2 * np) {
            // birth
            e_type = 1;
            e_pop = event_idx - np;
        } else if(event_idx < 3*np) {
            // death
            e_type = 2;
            e_pop = event_idx - 2*np;
        } else if(event_idx < 4*np) {
            e_type = 3;  // speciation
            e_pop = event_idx - 3*np;
        } else {
            e_type = 4;  // immigration
        }

        // record this simulation
        output(i, 0) = event_idx; // event ID
        output(i, 1) = e_type;    // event type
        output(i, 2) = e_pop;     // event pop
    }

    return output;
}
