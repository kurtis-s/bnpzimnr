#ifndef DPM_H
#define DPM_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fnStickBreak(NumericVector weights) {
    const int L = weights.size() + 1;
    NumericVector psi(L);
    NumericVector cum_prod(L);

    psi[0] = weights[0];
    cum_prod[0] = 1.0;
    for(int l=1; l<(L-1); l++) {
        cum_prod[l] = cum_prod[l-1]*(1.0-weights[l-1]);
        psi[l] = weights[l] * cum_prod[l];
    }
    psi[L-1] = std::max(1.0-sum(psi), 0.0); // So that we don't get negative probabilities (numerical stability)

    return psi;
}

// [[Rcpp::export]]
NumericMatrix fnSampStarMatCommon(NumericMatrix Param_mat, NumericMatrix C_mat, NumericVector sigma_2_vec, const double & a_star, const double & b_star_2, const double & L_trunc, bool truncated, double lower, double upper) {
    const double prior_precision = 1.0/b_star_2;
    const int & K = Param_mat.nrow();
    const int & J = Param_mat.ncol();

    double precision_curr;
    double l_curr_idx1;
    NumericMatrix draw_precision(K, L_trunc);
    std::fill(draw_precision.begin(), draw_precision.end(), prior_precision); // From the prior
    NumericMatrix unscaled_draw_mean(K, L_trunc);
    std::fill(unscaled_draw_mean.begin(), unscaled_draw_mean.end(), prior_precision * a_star); // From the prior

    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            l_curr_idx1 = C_mat(k, j);
            if(!R_IsNA(l_curr_idx1)) {
                if(R_IsNA(Param_mat(k, j))) stop("NA parameter when cluster indicator not NA in fnSampStarMat");
                precision_curr = 1.0/sigma_2_vec[k];
                draw_precision(k, l_curr_idx1 - 1) += precision_curr;
                unscaled_draw_mean(k, l_curr_idx1 - 1) += precision_curr * Param_mat(k, j);
            }
        }
    }

    double draw_mean;
    double draw_sd;
    NumericMatrix Param_star_mat(K, L_trunc);
    for(int l=0; l<L_trunc; l++) {
        for(int k=0; k<K; k++) {
            draw_mean = unscaled_draw_mean(k, l)/draw_precision(k, l);
            draw_sd = sqrt(1.0/draw_precision(k, l));
            // if(truncated) {
            //     Param_star_mat(k, l) = RcppTN::rtn1(draw_mean, draw_sd, lower, upper);
            // }
            // else {
            //     Param_star_mat(k, l) = R::rnorm(draw_mean, draw_sd);
            // }
            if(truncated) stop("Truncation draw not implemented"); // If changing add RcppTN in LinkingTo and Depends in description, and add "// #include <RcppTN.h>" and "// [[Rcpp::depends(RcppTN)]]" to this file
            Param_star_mat(k, l) = R::rnorm(draw_mean, draw_sd);
        }
    }

    return Param_star_mat;
}

// [[Rcpp::export]]
NumericMatrix fnSampStarMat(NumericMatrix Param_mat, NumericMatrix C_mat, NumericVector sigma_2_vec, const double & a_star, const double & b_star_2, const double & L_trunc) {
    return fnSampStarMatCommon(Param_mat, C_mat, sigma_2_vec, a_star, b_star_2, L_trunc, false, R_NegInf, R_PosInf);
}

// [[Rcpp::export]]
NumericMatrix fnSampStarMatTrunc(NumericMatrix Param_mat, NumericMatrix C_mat, NumericVector sigma_2_vec, const double & a_star, const double & b_star_2, const double & L_trunc, double lower, double upper) {
    return fnSampStarMatCommon(Param_mat, C_mat, sigma_2_vec, a_star, b_star_2, L_trunc, true, lower, upper);
}

// [[Rcpp::export]]
NumericVector fnCProbs(const double & param_kj, NumericVector psi_vec, NumericVector Param_star_k_vec, const double & sigma_k_2) {
    const double & L = psi_vec.length();
    NumericVector log_probs(L);
    for(int l=0; l<L; l++) {
        log_probs[l] = log(psi_vec[l]) + R::dnorm(param_kj, Param_star_k_vec[l], sqrt(sigma_k_2), TRUE);
    }
    NumericVector probs = exp(log_probs - max(log_probs));

    return probs;
}

void fnSampC(const double & k, const double & j, NumericMatrix C_mat, NumericMatrix Param_mat, IntegerVector d_vec, NumericVector psi_vec, NumericMatrix Param_star_mat, NumericVector sigma_2_vec) {
    const int & L = Param_star_mat.ncol();

    int l_idx1_old;
    int l_idx1_new;
    NumericVector probs;
    if(!R_IsNA(Param_mat(k, j)) && !(Param_mat(k, j)==0)) { // High abundance case and sum(delta_kj) > 1
        probs = fnCProbs(Param_mat(k, j), psi_vec, Param_star_mat(k, _), sigma_2_vec[k]);
        l_idx1_old = C_mat(k, j);
        l_idx1_new = sample(L, 1, FALSE, probs)(0);

        // Record new cluster assignment
        C_mat(k, j) = l_idx1_new;

        // Modify the cluster counts in place
        d_vec[l_idx1_old-1] -= 1;
        d_vec[l_idx1_new-1] += 1;
    }
    else {
        C_mat(k, j) = NA_REAL; // Low abundance case delta_kj=0, or beta_kj=0 because sum(delta_kj)=1
    }
}

// [[Rcpp::export]]
NumericMatrix fnSampCMat(NumericMatrix C_mat, NumericMatrix Param_mat, IntegerVector d_vec, NumericVector psi_vec, NumericMatrix Param_star_mat, NumericVector sigma_2_vec) {
    const int & K = Param_mat.nrow();
    const int & J = Param_mat.ncol();

    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            fnSampC(k, j, C_mat, Param_mat, d_vec, psi_vec, Param_star_mat, sigma_2_vec);
        }
    }

    return C_mat;
}

#endif
