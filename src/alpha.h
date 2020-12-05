#ifndef ALPHA_H
#define ALPHA_H

#include <Rcpp.h>
#include "utilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
double fnAlphaijLogLik(const double & alpha_ij, NumericVector mu_ij_vec, const double & s_j,
                        const int & lambda_ij, const double & eta_ell,
                        const double & w_ell, const double & upsilon_alpha,
                        const double & u_2, NumericVector y_ij_vec, NumericVector delta_ij_vec) {
    double prior_part = 0.0;
    if(lambda_ij == 1) {
        prior_part = R::dnorm(alpha_ij, eta_ell, sqrt(u_2), TRUE);
    }
    else if(lambda_ij == 0) {
        prior_part = R::dnorm(alpha_ij, (upsilon_alpha - w_ell*eta_ell)/(1 - w_ell), sqrt(u_2), TRUE);
    }
    else {
        Rcpp::Rcout << "alpha_ij=" << alpha_ij << std::endl;
        Rcpp::Rcout << "lambda_ij=" << lambda_ij << std::endl;
        stop("Inconsistent value for lambda.ij in fnAlphaijLogLik");
    }

    double likelihood_part = fnRepLogLik(mu_ij_vec, s_j, y_ij_vec, delta_ij_vec);

    return prior_part + likelihood_part;
}

// [[Rcpp::export]]
NumericMatrix fnSampAlphaMat(NumericMatrix Alpha_mat, NumericMatrix Mu_mat,
                NumericVector s_vec, NumericVector lambda_vec, NumericVector c_vec,
                NumericVector eta_vec, NumericVector w_vec, NumericVector rep_K,
                const double & upsilon_alpha, const double & u_2, NumericMatrix Y_mat, NumericMatrix Delta_mat,
                const double & alpha_ij_proposal_sd) {
    const int & n = Alpha_mat.nrow();
    const int & J = Alpha_mat.ncol();

    int K_i=-1; // Number of replicates for that subject
    double curr_log_lik;
    double prop_log_lik;
    double alpha_ij_prop;
    int ell_idx0;
    int ij=0;
    int ik=0;
    int ik_start=0;
    for(int j=0; j<J; j++) {
        ik=0;
        for(int i=0; i<n; i++) {
            K_i = rep_K[i];

            if(!R_IsNA(Alpha_mat(i, j))) { // Only sample alpha.ij if it isn't NA
                NumericVector mu_ij_vec(K_i);
                NumericVector mu_ij_prop(K_i);
                NumericVector delta_ij_vec(K_i);
                NumericVector y_ij_vec(K_i);
                ik_start = ik;
                alpha_ij_prop = Alpha_mat(i, j) + R::rnorm(0.0, alpha_ij_proposal_sd);
                for(int k=0; k<K_i; k++) {
                    mu_ij_vec[k] = Mu_mat(ik, j);
                    // mu_ij_prop[k] = Mu_mat(ik, j)/exp(Alpha_mat(i, j))*exp(alpha_ij_prop);
                    mu_ij_prop[k] = exp(log(Mu_mat(ik, j)) - Alpha_mat(i, j) + alpha_ij_prop);
                    delta_ij_vec[k] = Delta_mat(ik, j);
                    y_ij_vec[k] = Y_mat(ik, j);
                    ik++;
                }

                ell_idx0 = c_vec[ij] - 1;
                curr_log_lik = fnAlphaijLogLik(Alpha_mat(i, j), mu_ij_vec, s_vec[j], lambda_vec[ij], eta_vec[ell_idx0], w_vec[ell_idx0], upsilon_alpha, u_2, y_ij_vec, delta_ij_vec);
                prop_log_lik = fnAlphaijLogLik(alpha_ij_prop, mu_ij_prop, s_vec[j], lambda_vec[ij], eta_vec[ell_idx0], w_vec[ell_idx0], upsilon_alpha, u_2, y_ij_vec, delta_ij_vec);

                if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
                    Alpha_mat(i, j) = alpha_ij_prop;
                    for(int k=0; k<K_i; k++) Mu_mat(ik_start + k, j) = mu_ij_prop[k];
                }
            }
            else {
                // ik += K; // Still need to increment ik even if alpha.ij is NA
                for(int k=0; k<K_i; k++) ik++;
            }

            ij++;
        }
    }

    return Alpha_mat;
}

#endif
