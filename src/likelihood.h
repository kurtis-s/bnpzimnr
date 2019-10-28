#ifndef LIKELIHOOD_H
#define LIKELIHOOD_H

#include <Rcpp.h>
using namespace Rcpp;

/* Likelihood ############################################################### */
// [[Rcpp::export]]
double fnLogLikNoGamma(const double & y_ikj, const double & s_j, const double & mu_ikj, const double & delta_ikj) {
    // Does not include the y_ikj factorial part or the gamma part s_j
    const double s_j_inv = 1/s_j;
    const double param_product = mu_ikj * s_j;

    double ret;
    if(delta_ikj == 0) { // Low abundance case
        if(y_ikj != 0) stop("Inconsistent delta_ikj value");
        ret = 0.0;
    }
    else {
        // For numerical stability when param_product is infinite
        if(y_ikj != 0) { // High abundance case
            ret = -y_ikj * log(1.0/param_product + 1) - s_j_inv * log(1.0 + param_product);
        }
        else { // The above fails when y_ikj is 0 and param_product is infinite: 0*log(Inf) = Nan
            ret = -s_j_inv * log(1.0 + param_product);
        }
    }

    return ret;
}

// [[Rcpp::export]]
double fnLogLikNoFac(const double & y_ikj, const double & s_j,
                     const double & mu_ikj, const double & delta_ikj) {
    // Does not include the y_ikj factorial part
    const double s_j_inv = 1.0/s_j;

    return lgamma(y_ikj + s_j_inv) - lgamma(s_j_inv) +
        fnLogLikNoGamma(y_ikj, s_j, mu_ikj, delta_ikj);
}

// [[Rcpp::export]]
double fnOtuLogLikNoFac(const double & s_j, NumericVector mu_j_vec,
                        NumericVector rep_K, NumericVector y_j_vec, NumericVector delta_j_vec) {
    // Log-likelihood for the jth OTU over subject and replicates
    const int & n = rep_K.size();

    double total = 0;
    int ik=0;
    for(int i=0; i<n; i++) {
        for(int k=0; k<rep_K(i); k++) {
            total += fnLogLikNoFac(y_j_vec[ik], s_j, mu_j_vec[ik], delta_j_vec[ik]);
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
        log_lik_total += fnLogLikNoGamma(y_ik_vec[j], s_vec[j], mu_ik_vec[j], delta_ik_vec[j]);
    }

    return log_lik_total;
}

// [[Rcpp::export]]
double fnRepLogLikNoGamma(NumericVector mu_ij_vec, const double & s_j, NumericVector y_ij_vec, NumericVector delta_ij_vec) {
    // Log-likelihood for ij over replicates
    const int & K = delta_ij_vec.size();
    // For safety, can be deleted
    if(!(mu_ij_vec.size() == K && y_ij_vec.size() == K)) stop("Unexpected dimensions for fnRepLogLikNoGamma");

    double log_lik_total = 0.0;
    for(int k=0; k<K; k++) {
        log_lik_total += fnLogLikNoGamma(y_ij_vec[k], s_j, mu_ij_vec(k), delta_ij_vec[k]);
    }

    return log_lik_total;
}

// [[Rcpp::export]]
double fnSubjectLogLikNoGamma(const double & k, NumericVector y_j_vec, const double & s_j, NumericVector mu_j_vec, NumericVector delta_j_vec, NumericVector condition) {
    // Likelihood for kj over the subject i
    const int & N = y_j_vec.size();

    double lik = 0.0;
    for(int ik=0; ik<N; ik++) {
        if(condition[ik] == (k+1)) {
            lik += fnLogLikNoGamma(y_j_vec[ik], s_j, mu_j_vec[ik], delta_j_vec[ik]);
        }
    }

    return lik;
}

#endif
