#ifndef BETA_H
#define BETA_H

#include <Rcpp.h>
#include "likelihood.h"
#include "utilities.h"
using namespace Rcpp;

double fnXTBeta(NumericVector x_i, NumericVector beta_vec) {
    if ( x_i.size() != beta_vec.size() ) Rcpp::stop("x_i and beta_j dimension mismatch");

    double x_t_beta = 0.0;
    for(int k=0; k<x_i.size(); k++) {
        x_t_beta += x_i[k] * beta_vec[k];
    }

    return x_t_beta;
}

NumericVector fnGetMuJVecBetaProp(NumericMatrix X_mat, NumericVector mu_j_vec_curr, NumericVector delta_j_vec, NumericVector beta_j_vec_curr, NumericVector beta_j_vec_prop) {
    const int & N = X_mat.nrow();
    NumericVector mu_j_vec_prop(N);
    for(int i=0; i<N; i++) {
        if(delta_j_vec[i] == 1) {
            mu_j_vec_prop[i] = exp(log(mu_j_vec_curr[i]) - fnXTBeta(X_mat(i, _), beta_j_vec_curr) + fnXTBeta(X_mat(i, _), beta_j_vec_prop));
        } else if(delta_j_vec[i] == 0) {
            mu_j_vec_prop[i] = 0.0;
        }
        else {
            Rcpp::stop("Unexpected delta_i value");
        }
    }

    return mu_j_vec_prop;
}

NumericVector fnGetBetaJProp(NumericVector beta_j_vec, NumericVector beta_j_proposal_sd_vec) {
    const int & P = beta_j_vec.size();

    NumericVector beta_j_vec_prop(P);
    for(int p=0; p<P; p++) {
        beta_j_vec_prop[p] = beta_j_vec[p] + R::rnorm(0.0, beta_j_proposal_sd_vec[p]);
    }

    return beta_j_vec_prop;
}

double fnBetaPriorDens(NumericVector beta_j_vec, NumericVector tau_2_vec) {
    const int & P = beta_j_vec.size();
    double prior_part = 0.0;
    for(int p=0; p<P; p++) {
        prior_part += R::dnorm(beta_j_vec[p], 0.0, sqrt(tau_2_vec[p]), TRUE);
    }

    return prior_part;
}

// [[Rcpp::export]]
NumericMatrix fnSampBetaMat(NumericMatrix Beta_mat, NumericMatrix Mu_mat, NumericMatrix Delta_mat, NumericVector s_vec, NumericVector tau_2_vec, NumericMatrix X_mat, NumericMatrix Y_mat, NumericVector rep_K, NumericVector beta_j_proposal_sd_vec) {
    const int & J = Beta_mat.ncol();
    NumericVector mu_j_vec_prop;
    NumericVector beta_j_vec_prop;
    double prior_part_curr;
    double prior_part_prop;
    double lik_part_curr;
    double lik_part_prop;
    double curr_log_lik;
    double prop_log_lik;

    for(int j=0; j<J; j++) {
        beta_j_vec_prop = fnGetBetaJProp(Beta_mat(_, j), beta_j_proposal_sd_vec);
        mu_j_vec_prop = fnGetMuJVecBetaProp(X_mat, Mu_mat(_, j), Delta_mat(_, j), Beta_mat(_, j), beta_j_vec_prop);

        prior_part_curr = fnBetaPriorDens(Beta_mat(_, j), tau_2_vec);
        prior_part_prop = fnBetaPriorDens(beta_j_vec_prop, tau_2_vec);
        lik_part_curr = fnOtuLogLik(s_vec[j], Mu_mat(_, j), rep_K, Y_mat(_, j), Delta_mat(_, j));
        lik_part_prop = fnOtuLogLik(s_vec[j], mu_j_vec_prop, rep_K, Y_mat(_, j), Delta_mat(_, j));
        curr_log_lik = prior_part_curr + lik_part_curr;
        prop_log_lik = prior_part_prop + lik_part_prop;

        if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
            Beta_mat(_, j) = beta_j_vec_prop;
            Mu_mat(_, j) = mu_j_vec_prop;
        }
    }

    return Beta_mat;
}

#endif
