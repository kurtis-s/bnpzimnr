#ifndef BETA_H
#define BETA_H

#include <Rcpp.h>
#include "likelihood.h"
#include "utilities.h"
#include "dpm.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fnGetMuJVecProp(const double & k, NumericVector mu_j_vec_curr, const double & beta_kj_curr, const double & beta_kj_prop, NumericVector r_vec, NumericVector alpha_j_vec, NumericVector condition, NumericVector subject) {
    const double & N = mu_j_vec_curr.size();

    int i_idx0;
    NumericVector mu_j_vec_prop(N);
    for(int ik=0; ik<N; ik++) {
        if(condition[ik] == (k+1)) {
            i_idx0 = subject[ik] - 1;
            if(R_IsNA(alpha_j_vec[i_idx0])) {
                mu_j_vec_prop[ik] = 0.0;
            }
            else if(R_IsNA(beta_kj_curr) && !R_IsNA(beta_kj_prop)) { // In this case we must calculate mu.ikj
                mu_j_vec_prop[ik] = exp(r_vec[ik] + alpha_j_vec[i_idx0] + beta_kj_prop);
            }
            else if(R_IsNA(beta_kj_prop)) {
                mu_j_vec_prop[ik] = 0.0;
            }
            else {
                mu_j_vec_prop[ik] = exp(log(mu_j_vec_curr[ik]) - beta_kj_curr + beta_kj_prop);
            }
        }
        else {
            mu_j_vec_prop[ik] = mu_j_vec_curr[ik];
        }
    }

    return mu_j_vec_prop;
}

void fnSampBetakj(const double & k, const double & j, const double & beta_kj_curr, const double & c_kj, NumericVector delta_j_vec, const double & sigma_2_beta_k, NumericMatrix Beta_star_mat, NumericVector r_vec, NumericVector alpha_j_vec, NumericMatrix Beta_mat, NumericMatrix Mu_mat, const double & s_j, NumericVector y_j_vec, NumericVector condition, NumericVector subject, const double & beta_kj_proposal_sd) {
    double beta_kj_prop = Beta_mat(k, j) + R::rnorm(0.0, beta_kj_proposal_sd);
    double prior_part_curr = R::dnorm(beta_kj_curr, Beta_star_mat(k, c_kj-1), sqrt(sigma_2_beta_k), TRUE);
    double prior_part_prop = R::dnorm(beta_kj_prop, Beta_star_mat(k, c_kj-1), sqrt(sigma_2_beta_k), TRUE);

    double lik_part_curr = fnSubjectLogLikNoGamma(k, y_j_vec, s_j, Mu_mat(_, j), delta_j_vec, condition);
    NumericVector mu_j_vec_prop = fnGetMuJVecProp(k, Mu_mat(_, j), beta_kj_curr, beta_kj_prop, r_vec, alpha_j_vec, condition, subject);

    double lik_part_prop = fnSubjectLogLikNoGamma(k, y_j_vec, s_j, mu_j_vec_prop, delta_j_vec, condition);

    double curr_log_lik = prior_part_curr + lik_part_curr;
    double prop_log_lik = prior_part_prop + lik_part_prop;

    if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
        Beta_mat(k, j) = beta_kj_prop;
        Mu_mat(_, j) = mu_j_vec_prop;
    }
}

// [[Rcpp::export]]
NumericMatrix fnSampBetaMat(NumericMatrix Beta_mat, NumericMatrix Mu_mat, NumericMatrix C_beta_mat, NumericMatrix Beta_star_mat, NumericMatrix Delta_mat, NumericVector sigma_2_beta_vec, NumericVector s_vec, NumericVector r_vec, NumericMatrix Alpha_mat, NumericMatrix Y_mat, NumericVector condition, NumericVector subject, NumericMatrix B_mat, NumericMatrix Delta_sum_mat, const double & beta_kj_proposal_sd) {
    const int & K = Beta_mat.nrow(); // Only works if the same number of categories for each subject
    const int & J = Beta_mat.ncol();

    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            if(B_mat(k, j)==1 && Delta_sum_mat(k, j)>=1) {
                fnSampBetakj(k, j, Beta_mat(k, j), C_beta_mat(k, j), Delta_mat(_, j), sigma_2_beta_vec[k], Beta_star_mat, r_vec, Alpha_mat(_, j), Beta_mat, Mu_mat, s_vec[j], Y_mat(_, j), condition, subject, beta_kj_proposal_sd);
            }
        }
    }

    return(Beta_mat);
}

// [[Rcpp::export]]
NumericMatrix fnSampBetaMatMarg(NumericMatrix Beta_mat, NumericMatrix Mu_mat, NumericVector psi_beta_vec, NumericMatrix Beta_star_mat, NumericMatrix Delta_mat, NumericVector sigma_2_beta_vec, NumericVector s_vec, NumericVector r_vec, NumericMatrix Alpha_mat, NumericMatrix Y_mat, NumericVector condition, NumericVector subject, NumericMatrix B_mat, NumericMatrix Delta_sum_mat, const double & beta_kj_proposal_sd) {
    const int & K = Beta_mat.nrow(); // Only works if the same number of categories for each subject
    const int & J = Beta_mat.ncol();
    const int & L = psi_beta_vec.length();

    double beta_kj_prop;
    double prior_part_curr;
    double prior_part_prop;
    double lik_part_curr;
    double lik_part_prop;
    double curr_log_lik;
    double prop_log_lik;
    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            if(B_mat(k, j)==1 && Delta_sum_mat(k, j)>=1) {
                beta_kj_prop = Beta_mat(k, j) + R::rnorm(0.0, beta_kj_proposal_sd);

                prior_part_curr = 0.0;
                prior_part_prop = 0.0;
                for(int l=0; l<L; l++) {
                    prior_part_curr += psi_beta_vec[l] * R::dnorm(Beta_mat(k, j), Beta_star_mat(k, l), sqrt(sigma_2_beta_vec[k]), FALSE);
                    prior_part_prop += psi_beta_vec[l] * R::dnorm(beta_kj_prop, Beta_star_mat(k, l), sqrt(sigma_2_beta_vec[k]), FALSE);
                }
                prior_part_curr = log(prior_part_curr);
                prior_part_prop = log(prior_part_prop);

                lik_part_curr = fnSubjectLogLikNoGamma(k, Y_mat(_, j), s_vec[j], Mu_mat(_, j), Delta_mat(_, j), condition);
                NumericVector mu_j_vec_prop = fnGetMuJVecProp(k, Mu_mat(_, j), Beta_mat(k, j), beta_kj_prop, r_vec, Alpha_mat(_, j), condition, subject);

                lik_part_prop = fnSubjectLogLikNoGamma(k, Y_mat(_, j), s_vec[j], mu_j_vec_prop, Delta_mat(_, j), condition);

                curr_log_lik = prior_part_curr + lik_part_curr;
                prop_log_lik = prior_part_prop + lik_part_prop;

                if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
                    Beta_mat(k, j) = beta_kj_prop;
                    Mu_mat(_, j) = mu_j_vec_prop;
                }
            }
        }
    }

    return(Beta_mat);
}


#endif
