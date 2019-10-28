#ifndef NORMALIZATION_H
#define NORMALIZATION_H

#include <Rcpp.h>
#include "utilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
double fnRtkLogLik(const double & r_ik, NumericVector mu_ik_vec, NumericVector s_vec,
                        const int & lambda_ik, const double & eta_ell,
                        const double & w_ell, const double & upsilon_r,
                        const double & u_2, NumericVector y_ik_vec, NumericVector delta_ik_vec) {
    double prior_part = 0.0;
    if(lambda_ik == 1) {
        prior_part = R::dnorm(r_ik, eta_ell, sqrt(u_2), TRUE);
    }
    else { // lambda_ik == 0
        prior_part = R::dnorm(r_ik, (upsilon_r - w_ell*eta_ell)/(1 - w_ell), sqrt(u_2), TRUE);
    }

    double likelihood_part = fnSampleLogLikNoGamma(mu_ik_vec, s_vec, y_ik_vec, delta_ik_vec);

    return prior_part + likelihood_part;
}

// [[Rcpp::export]]
NumericVector fnSampRVec(NumericVector r_vec, NumericMatrix Mu_mat,
                NumericVector s_vec, NumericVector lambda_vec, NumericVector c_vec,
                NumericVector eta_vec, NumericVector w_vec, NumericVector rep_K,
                const double & upsilon_r, const double & u_2, NumericMatrix Y_mat, NumericMatrix Delta_mat,
                NumericVector r_vec_proposal_sd) {
    const int & N = r_vec.size();

    int ell_idx0;
    double curr_log_lik;
    double prop_log_lik;
    double r_ik_prop;
    NumericVector mu_ik_prop;
    for(int ik=0; ik<N; ik++) {
        ell_idx0 = c_vec[ik] - 1;

        r_ik_prop = r_vec[ik] + R::rnorm(0.0, r_vec_proposal_sd[ik]);
        mu_ik_prop = exp(log(Mu_mat(ik, _)) - r_vec[ik] + r_ik_prop);

        curr_log_lik = fnRtkLogLik(r_vec[ik], Mu_mat(ik, _), s_vec, lambda_vec[ik], eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(ik, _), Delta_mat(ik, _));
        prop_log_lik = fnRtkLogLik(r_ik_prop, mu_ik_prop, s_vec, lambda_vec[ik], eta_vec(ell_idx0), w_vec(ell_idx0), upsilon_r, u_2, Y_mat(ik, _), Delta_mat(ik, _));

        if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
            r_vec[ik] = r_ik_prop;
            Mu_mat(ik, _) = mu_ik_prop;
        }
    }

    return r_vec;
}

#endif
