#ifndef EPSILON_H
#define EPSILON_H

#include <Rcpp.h>
#include "utilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix fnSampXiMat(NumericMatrix Xi_mat, NumericMatrix Epsilon_mat, NumericMatrix Delta_mat, NumericMatrix Xi_star_mat, NumericMatrix C_xi_mat, NumericVector sigma_2_xi_vec, NumericVector condition, const double & xi_kj_proposal_sd) {
    const double & K= Xi_mat.nrow();
    const double & J = Xi_mat.ncol();
    const double & N = Delta_mat.nrow();

    double prior_part_curr;
    double prior_part_prop;
    double lik_part_curr;
    double lik_part_prop;
    double curr_log_lik;
    double prop_log_lik;
    double xi_kj_prop;
    double epsilon_kj_prop;
    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            xi_kj_prop = Xi_mat(k, j) + R::rnorm(0, xi_kj_proposal_sd);
            epsilon_kj_prop = R::pnorm(xi_kj_prop, 0.0, 1.0, true, false);

            prior_part_curr = R::dnorm(Xi_mat(k, j), Xi_star_mat(k, C_xi_mat(k, j)-1.0), sqrt(sigma_2_xi_vec[k]), TRUE);
            prior_part_prop = R::dnorm(xi_kj_prop, Xi_star_mat(k, C_xi_mat(k, j)-1.0), sqrt(sigma_2_xi_vec[k]), TRUE);
            lik_part_curr = 0.0;
            lik_part_prop = 0.0;
            for(int ik=0; ik<N; ik++) {
                if(k==(condition[ik]-1)) {
                    lik_part_curr += (1-Delta_mat(ik, j)) * log(Epsilon_mat(k, j)) + Delta_mat(ik, j) * log1p(-Epsilon_mat(k, j));
                    lik_part_prop += (1-Delta_mat(ik, j)) * log(epsilon_kj_prop) + Delta_mat(ik, j) * log1p(-epsilon_kj_prop);
                }
            }

            curr_log_lik = prior_part_curr + lik_part_curr;
            prop_log_lik = prior_part_prop + lik_part_prop;

            if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
                Xi_mat(k, j) = xi_kj_prop;
                Epsilon_mat(k, j) = epsilon_kj_prop;
            }
        }
    }

    return Xi_mat;
}

// [[Rcpp::export]]
NumericMatrix fnSampXiMatMarg(NumericMatrix Xi_mat, NumericMatrix Epsilon_mat, NumericMatrix Delta_mat, NumericMatrix Xi_star_mat, NumericVector psi_xi_vec, NumericVector sigma_2_xi_vec, NumericVector condition, const double & xi_kj_proposal_sd) {
    const double & K= Xi_mat.nrow();
    const double & J = Xi_mat.ncol();
    const double & N = Delta_mat.nrow();
    const int & L = psi_xi_vec.length();

    double prior_part_curr;
    double prior_part_prop;
    double lik_part_curr;
    double lik_part_prop;
    double curr_log_lik;
    double prop_log_lik;
    double xi_kj_prop;
    double epsilon_kj_prop;
    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            xi_kj_prop = Xi_mat(k, j) + R::rnorm(0, xi_kj_proposal_sd);
            epsilon_kj_prop = R::pnorm(xi_kj_prop, 0.0, 1.0, true, false);

            prior_part_curr = 0.0;
            prior_part_prop = 0.0;
            for(int l=0; l<L; l++) {
                prior_part_curr += psi_xi_vec[l] * R::dnorm(Xi_mat(k, j), Xi_star_mat(k, l), sqrt(sigma_2_xi_vec[k]), FALSE);
                prior_part_prop += psi_xi_vec[l] * R::dnorm(xi_kj_prop, Xi_star_mat(k, l), sqrt(sigma_2_xi_vec[k]), FALSE);
            }
            prior_part_curr = log(prior_part_curr);
            prior_part_prop = log(prior_part_prop);


            lik_part_curr = 0.0;
            lik_part_prop = 0.0;
            for(int ik=0; ik<N; ik++) {
                if(k==(condition[ik]-1)) {
                    lik_part_curr += (1-Delta_mat(ik, j)) * log(Epsilon_mat(k, j)) + Delta_mat(ik, j) * log1p(-Epsilon_mat(k, j));
                    lik_part_prop += (1-Delta_mat(ik, j)) * log(epsilon_kj_prop) + Delta_mat(ik, j) * log1p(-epsilon_kj_prop);
                }
            }

            curr_log_lik = prior_part_curr + lik_part_curr;
            prop_log_lik = prior_part_prop + lik_part_prop;

            if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
                Xi_mat(k, j) = xi_kj_prop;
                Epsilon_mat(k, j) = epsilon_kj_prop;
            }
        }
    }

    return Xi_mat;
}

#endif
