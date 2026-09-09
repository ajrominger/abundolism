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
                          int nstep) {

    vec m = m_prop % g;
    int nsim = la.n_elem;

    Rcpp::List out(nsim);

    int total_events = 4 * np + 1;

    uvec x(np);
    vec s(np);
    vec stau(np);
    uvec sfull(np);
    double tt;
    bool any_full_spec = false;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int> sample_pop(0, np - 1);

    for(int r = 0; r < nsim; r++) {

        x.zeros();
        s.zeros();
        stau.zeros();
        sfull.zeros();
        tt = 0.0;
        any_full_spec = false;

        // pre-allocate this sim's matrix at the worst case (nstep rows),
        // trim to actual size before storing in the list
        mat sim_output(nstep, 4);
        int step_count = 0;

        for(int i = 0; i < nstep; i++) {

            vec this_la = conv_to<vec>::from(x) * la(r);
            vec this_mu = conv_to<vec>::from(x) * mu(r);
            vec this_m = conv_to<vec>::from(x) * m(r);
            vec this_nu = conv_to<vec>::from(x) * nu(r);
            double this_g = g(r);

            std::vector<double> probs(total_events);
            for(int j = 0; j < np; j++) {
                probs[j] = this_m(j);
                probs[np + j] = this_la(j);
                probs[2*np + j] = this_mu(j);
                probs[3*np + j] = this_nu(j);
            }
            probs[4*np] = this_g;

            std::discrete_distribution<int> sample_event(probs.begin(), probs.end());
            int event_idx = sample_event(gen);

            int e_type;
            int e_pop = -1;

            if(event_idx < np) {
                e_type = 0;
                e_pop = event_idx;
            } else if(event_idx < 2 * np) {
                e_type = 1;
                e_pop = event_idx - np;
            } else if(event_idx < 3 * np) {
                e_type = 2;
                e_pop = event_idx - 2 * np;
            } else if(event_idx < 4 * np) {
                e_type = 3;
                e_pop = event_idx - 3 * np;
            } else {
                e_type = 4;
            }

            double total_rate = std::accumulate(probs.begin(), probs.end(), 0.0);
            std::exponential_distribution<double> exp_dist(total_rate);
            double st = exp_dist(gen);

            for(int j = 0; j < np; j++) {
                if(s(j) > 0) {
                    s(j) += st;
                }
            }

            if(e_type == 0 || e_type == 4) {
                int receiving_pop = sample_pop(gen);
                x(receiving_pop) ++;
                if(s(receiving_pop) > 0) {
                    stau(receiving_pop) += xi(r) / x(receiving_pop);
                }
            } else if(e_type == 1) {
                x(e_pop) ++;
            } else if(e_type == 2) {
                x(e_pop) --;
            } else if(e_type == 3) {
                if(s(e_pop) == 0) {
                    s(e_pop) = st;
                }
            }

            tt += st;

            for(int j = 0; j < np; j++) {
                if(x(j) == 0) {
                    stau(j) = tau(r);
                    s(j) = 0;
                }
            }

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
