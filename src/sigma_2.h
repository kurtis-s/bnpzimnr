#ifndef SIGMA_2_H
#define SIGMA_2_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector fnSampSigma2(NumericMatrix Param_mat, NumericMatrix Param_star_mat, NumericMatrix C_mat, const double & a_param, const double & b_param) {
    const int & K = Param_mat.nrow();
    const int & J = Param_mat.ncol();


    NumericVector a_prime(K);
    std::fill(a_prime.begin(), a_prime.end(), a_param);
    NumericVector b_prime(K);
    std::fill(b_prime.begin(), b_prime.end(), b_param);
    for(int j=0; j<J; j++) {
        for(int k=0; k<K; k++) {
            if( (Param_mat(k, j) != 0.0) && (!R_IsNA(Param_mat(k, j))) ) { // To handle the case where theta_kj doesn't exist
                a_prime[k] += 1.0/2.0;
                b_prime[k] += (1.0/2.0) * pow(Param_mat(k, j) - Param_star_mat(k, C_mat(k, j)-1), 2.0);
            }
        }
    }

    NumericVector sigma_vec(K);
    for(int k=0; k<K; k++) {
        sigma_vec[k] = (1.0/rgamma(1.0, a_prime[k], 1.0/b_prime[k]))[0];
    }

    return sigma_vec;
}

#endif
