#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <Rcpp.h>
using namespace Rcpp;

/* Likelihood ############################################################### */
// [[Rcpp::export]]
double fnLogLik(const double & y_ikj, const double & s_j, const double & mu_ikj, const double & delta_ikj) {
    double ret;
    if(delta_ikj == 0) {
        ret = 0.0;
    } else {
        ret = R::dnbinom_mu(y_ikj, 1.0/s_j, mu_ikj, TRUE);
    }

    return ret;
}

// [[Rcpp::export]]
double fnOtuLogLik(const double & s_j, NumericVector mu_j_vec,
                        NumericVector rep_K, NumericVector y_j_vec, NumericVector delta_j_vec) {
    // Log-likelihood for the jth OTU over subject and replicates
    const int & n = rep_K.size();

    double total = 0.0;
    int ik=0;
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            total += fnLogLik(y_j_vec[ik], s_j, mu_j_vec[ik], delta_j_vec[ik]);
            ik++;
        }
    }

    return total;
}

// [[Rcpp::export]]
double fnSampleLogLikNoGamma(NumericVector mu_ik_vec, NumericVector s_vec,
                             NumericVector y_ik_vec, NumericVector delta_ik_vec) {
    // Log-likelihood for sample ik over OTUs
    const int & J = mu_ik_vec.size();

    double log_lik_total = 0.0;
    for(int j=0; j<J; j++) {
        log_lik_total += fnLogLik(y_ik_vec[j], s_vec[j], mu_ik_vec[j], delta_ik_vec[j]);
    }

    return log_lik_total;
}

// [[Rcpp::export]]
double fnRepLogLik(NumericVector mu_ij_vec, const double & s_j, NumericVector y_ij_vec, NumericVector delta_ij_vec) {
    // Log-likelihood for ij over replicates
    const int & K_i = delta_ij_vec.size();
    // For safety, can be deleted
    if(!(mu_ij_vec.size() == K_i && y_ij_vec.size() == K_i)) stop("Unexpected dimensions for fnRepLogLik");

    double log_lik_total = 0.0;
    for(int k=0; k<K_i; k++) {
        log_lik_total += fnLogLik(y_ij_vec[k], s_j, mu_ij_vec(k), delta_ij_vec[k]);
    }

    return log_lik_total;
}

// [[Rcpp::export]]
double fnSubjectLogLik(const double & k, NumericVector y_j_vec, const double & s_j, NumericVector mu_j_vec, NumericVector delta_j_vec, NumericVector condition) {
    // Likelihood for kj over the subject i
    const int & N = y_j_vec.size();

    double lik = 0.0;
    for(int ik=0; ik<N; ik++) {
        if(condition[ik] == (k+1)) {
            lik += fnLogLik(y_j_vec[ik], s_j, mu_j_vec[ik], delta_j_vec[ik]);
        }
    }

    return lik;
}

#endif
