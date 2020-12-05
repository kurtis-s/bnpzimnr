#ifndef THETA_H
#define THETA_H

#include <Rcpp.h>
#include "likelihood.h"
#include "utilities.h"
#include "dpm.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fnGetMuJVecProp(const double & k, NumericVector mu_j_vec_curr, const double & theta_kj_curr, const double & theta_kj_prop, NumericVector r_vec, NumericVector alpha_j_vec, NumericVector condition, NumericVector subject) {
    const double & N = mu_j_vec_curr.size();

    int i_idx0;
    NumericVector mu_j_vec_prop(N);
    for(int ik=0; ik<N; ik++) {
        if(condition[ik] == (k+1)) {
            i_idx0 = subject[ik] - 1;
            if(R_IsNA(alpha_j_vec[i_idx0])) {
                mu_j_vec_prop[ik] = 0.0;
            }
            else if(R_IsNA(theta_kj_curr) && !R_IsNA(theta_kj_prop)) { // In this case we must calculate mu.ikj
                mu_j_vec_prop[ik] = exp(r_vec[ik] + alpha_j_vec[i_idx0] + theta_kj_prop);
            }
            else if(R_IsNA(theta_kj_prop)) {
                mu_j_vec_prop[ik] = 0.0;
            }
            else {
                // mu_j_vec_prop[ik] = mu_j_vec_curr[ik]/exp(theta_kj_curr)*exp(theta_kj_prop);
                mu_j_vec_prop[ik] = exp(log(mu_j_vec_curr[ik]) - theta_kj_curr + theta_kj_prop);
            }
        }
        else {
            mu_j_vec_prop[ik] = mu_j_vec_curr[ik];
        }
    }

    return mu_j_vec_prop;
}

void fnSampThetakj(const double & k, const double & j, const double & theta_kj_curr, const double & c_kj, NumericVector delta_j_vec, const double & sigma_2_theta_k, NumericMatrix Theta_star_mat, NumericVector r_vec, NumericVector alpha_j_vec, NumericMatrix Theta_mat, NumericMatrix Mu_mat, const double & s_j, NumericVector y_j_vec, NumericVector condition, NumericVector subject, const double & theta_kj_proposal_sd) {
    double theta_kj_prop = Theta_mat(k, j) + R::rnorm(0.0, theta_kj_proposal_sd);
    double prior_part_curr = R::dnorm(theta_kj_curr, Theta_star_mat(k, c_kj-1), sqrt(sigma_2_theta_k), TRUE);
    double prior_part_prop = R::dnorm(theta_kj_prop, Theta_star_mat(k, c_kj-1), sqrt(sigma_2_theta_k), TRUE);

    double lik_part_curr = fnSubjectLogLik(k, y_j_vec, s_j, Mu_mat(_, j), delta_j_vec, condition);
    NumericVector mu_j_vec_prop = fnGetMuJVecProp(k, Mu_mat(_, j), theta_kj_curr, theta_kj_prop, r_vec, alpha_j_vec, condition, subject);

    double lik_part_prop = fnSubjectLogLik(k, y_j_vec, s_j, mu_j_vec_prop, delta_j_vec, condition);

    double curr_log_lik = prior_part_curr + lik_part_curr;
    double prop_log_lik = prior_part_prop + lik_part_prop;

    if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
        Theta_mat(k, j) = theta_kj_prop;
        Mu_mat(_, j) = mu_j_vec_prop;
    }
}

// [[Rcpp::export]]
NumericMatrix fnSampThetaMat(NumericMatrix Theta_mat, NumericMatrix Mu_mat, NumericMatrix C_theta_mat, NumericMatrix Theta_star_mat, NumericMatrix Delta_mat, NumericVector sigma_2_theta_vec, NumericVector s_vec, NumericVector r_vec, NumericMatrix Alpha_mat, NumericMatrix Y_mat, NumericVector condition, NumericVector subject, NumericMatrix B_mat, NumericMatrix Delta_sum_mat, const double & theta_kj_proposal_sd) {
    const int & K = Theta_mat.nrow(); // Only works if the same number of categories for each subject
    const int & J = Theta_mat.ncol();

    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            if(B_mat(k, j)==1 && Delta_sum_mat(k, j)>=1) {
                fnSampThetakj(k, j, Theta_mat(k, j), C_theta_mat(k, j), Delta_mat(_, j), sigma_2_theta_vec[k], Theta_star_mat, r_vec, Alpha_mat(_, j), Theta_mat, Mu_mat, s_vec[j], Y_mat(_, j), condition, subject, theta_kj_proposal_sd);
            }
        }
    }

    return(Theta_mat);
}

// [[Rcpp::export]]
NumericMatrix fnSampThetaMatMarg(NumericMatrix Theta_mat, NumericMatrix Mu_mat, NumericVector psi_theta_vec, NumericMatrix Theta_star_mat, NumericMatrix Delta_mat, NumericVector sigma_2_theta_vec, NumericVector s_vec, NumericVector r_vec, NumericMatrix Alpha_mat, NumericMatrix Y_mat, NumericVector condition, NumericVector subject, NumericMatrix B_mat, NumericMatrix Delta_sum_mat, const double & theta_kj_proposal_sd) {
    const int & K = Theta_mat.nrow(); // Only works if the same number of categories for each subject
    const int & J = Theta_mat.ncol();
    const int & L = psi_theta_vec.length();

    double theta_kj_prop;
    double prior_part_curr;
    double prior_part_prop;
    double lik_part_curr;
    double lik_part_prop;
    double curr_log_lik;
    double prop_log_lik;
    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            if(B_mat(k, j)==1 && Delta_sum_mat(k, j)>=1) {
                theta_kj_prop = Theta_mat(k, j) + R::rnorm(0.0, theta_kj_proposal_sd);

                prior_part_curr = 0.0;
                prior_part_prop = 0.0;
                for(int l=0; l<L; l++) {
                    prior_part_curr += psi_theta_vec[l] * R::dnorm(Theta_mat(k, j), Theta_star_mat(k, l), sqrt(sigma_2_theta_vec[k]), FALSE);
                    prior_part_prop += psi_theta_vec[l] * R::dnorm(theta_kj_prop, Theta_star_mat(k, l), sqrt(sigma_2_theta_vec[k]), FALSE);
                }
                prior_part_curr = log(prior_part_curr);
                prior_part_prop = log(prior_part_prop);

                lik_part_curr = fnSubjectLogLik(k, Y_mat(_, j), s_vec[j], Mu_mat(_, j), Delta_mat(_, j), condition);
                NumericVector mu_j_vec_prop = fnGetMuJVecProp(k, Mu_mat(_, j), Theta_mat(k, j), theta_kj_prop, r_vec, Alpha_mat(_, j), condition, subject);

                lik_part_prop = fnSubjectLogLik(k, Y_mat(_, j), s_vec[j], mu_j_vec_prop, Delta_mat(_, j), condition);

                curr_log_lik = prior_part_curr + lik_part_curr;
                prop_log_lik = prior_part_prop + lik_part_prop;

                if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
                    Theta_mat(k, j) = theta_kj_prop;
                    Mu_mat(_, j) = mu_j_vec_prop;
                }
            }
        }
    }

    return(Theta_mat);
}


#endif
