#ifndef OVERDISPERSION_H
#define OVERDISPERSION_H

#include <Rcpp.h>
#include "utilities.h"
#include "likelihood.h"
using namespace Rcpp;

/* Overdispersion ########################################################### */
// ### s_j
// [[Rcpp::export]]
double fnSjTildeLogLik(const double & s_j, NumericVector mu_j_vec,
                       NumericVector rep_K, NumericVector y_j_vec, NumericVector delta_j_vec,
                       const double & h_scal, const double & kappa_2) {
    return R::dnorm(log(s_j), h_scal, sqrt(kappa_2), TRUE) +
        fnOtuLogLikNoFac(s_j, mu_j_vec, rep_K, y_j_vec, delta_j_vec);
}

// [[Rcpp::export]]
NumericVector fnSampSVec(NumericVector s_vec, NumericMatrix Mu_mat,
                NumericVector rep_K, NumericMatrix Y_mat, NumericMatrix Delta_mat,
                const double & h_scal, const double & kappa_2,
                NumericVector s_vec_proposal_sd) {
    const int & J = s_vec.size();

    double curr_log_lik = 0.0;
    double prop_log_lik = 0.0;
    double s_j_prop = 0.0;
    for(int j=0; j<J; j++) {
        s_j_prop = exp(log(s_vec(j)) + R::rnorm(0.0, s_vec_proposal_sd[j]));

        curr_log_lik = fnSjTildeLogLik(s_vec(j), Mu_mat(_, j), rep_K, Y_mat(_, j), Delta_mat(_, j), h_scal, kappa_2);
        prop_log_lik = fnSjTildeLogLik(s_j_prop, Mu_mat(_, j), rep_K, Y_mat(_, j), Delta_mat(_, j), h_scal, kappa_2);

        if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
            s_vec(j) = s_j_prop;
        }
    }

    return s_vec;
}

#endif
