#ifndef MEAN_CONSTRAINT_H
#define MEAN_CONSTRAINT_H

#include <Rcpp.h>
#include "utilities.h"
using namespace Rcpp;

// [[Rcpp::export]]
double fnWlLogLik(const int & ell_idx1, const double & w_ell,
                  NumericVector param_vec, NumericVector lambda_vec,
                  NumericVector c_vec, NumericVector eta_vec,
                  const double & upsilon, const double & a_w,
                  const double & b_w, const double & u_2) {
    // Guard to make sure w_ell is in the support
    if(w_ell < 0.0 || w_ell > 1.0) {
        return R_NegInf;
    }

    double normal_mean = (upsilon - w_ell*eta_vec[ell_idx1-1])/(1.0 - w_ell);

    double w_part_total = 0.0; // The part of the full conditional from w
    double param_part_total = 0.0; // The part of the full conditional from r_ik or theta0_j
    double lambda_cur;
    double c_cur;
    double param_cur;
    for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is tk for w_\ell^r; or j for w_\ell^\theta
        param_cur = param_vec[loop_idx];

        if(!R_IsNA(param_cur)) {
            lambda_cur = lambda_vec[loop_idx];
            c_cur = c_vec[loop_idx];

            if(c_cur == ell_idx1) {
                if(lambda_cur == 0) {
                    w_part_total += log(1.0 - w_ell);
                    param_part_total += R::dnorm(param_cur, normal_mean, sqrt(u_2), TRUE);
                }
                else { // lambda = 1
                    w_part_total += log(w_ell);
                }
            }
        }
    }

    double prior_part = (a_w - 1.0)*log(w_ell) + (b_w - 1.0)*log(1.0 - w_ell); // From the prior

    return prior_part + w_part_total + param_part_total;
}

// [[Rcpp::export]]
NumericVector fnSampWVec(NumericVector w_vec, NumericVector param_vec,
                NumericVector lambda_vec, NumericVector c_vec,
                NumericVector eta_vec, const double & upsilon,
                const double & a_w, const double & b_w, const double & u_2,
                NumericVector w_vec_proposal_sd) {
    const int L = w_vec.size();

    double jacob_curr;
    double jacob_prop;
    double w_ell_curr_logit;
    double w_ell_prop_logit;
    double curr_log_lik = 0.0;
    double prop_log_lik = 0.0;
    double w_ell_prop = -1.0;
    for(int ell=0; ell<L; ell++) {
        w_ell_curr_logit = log(w_vec[ell]/(1.0-w_vec[ell]));
        w_ell_prop_logit = w_ell_curr_logit + R::rnorm(0.0, w_vec_proposal_sd[ell]);
        w_ell_prop = 1.0/(1.0+exp(w_ell_prop_logit));

        jacob_curr = -w_ell_curr_logit - 2.0*log(exp(-w_ell_curr_logit) + 1.0);
        jacob_prop = -w_ell_prop_logit - 2.0*log(exp(-w_ell_prop_logit) + 1.0);

        curr_log_lik = fnWlLogLik(ell + 1, w_vec(ell), param_vec, lambda_vec, c_vec, eta_vec, upsilon, a_w, b_w, u_2) + jacob_curr;
        prop_log_lik = fnWlLogLik(ell + 1, w_ell_prop, param_vec, lambda_vec, c_vec, eta_vec, upsilon, a_w, b_w, u_2) + jacob_prop;

        if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
            w_vec(ell) = w_ell_prop;
        }
    }

    return w_vec;
}

// [[Rcpp::export]]
NumericVector fnSampEtaVec(NumericVector eta_vec, NumericVector param_vec, NumericVector w_vec,
                  NumericVector lambda_vec, NumericVector c_vec,
                  const double & u_2, const double & upsilon, const double & b_eta_2) {
    const int L = eta_vec.size();
    const double & m3 = upsilon;
    const double rho3 = 1/b_eta_2;

    int ell_idx1;
    double norm_mean;
    double norm_precision;
    double norm_variance;
    double rho1;
    double rho2;
    double m1;
    double m2;
    double m1_total = 0.0;
    double m2_total = 0.0;
    double set_A_size = 0.0;
    double set_B_size = 0.0;
    for(int ell=0; ell<L; ell++) {
        ell_idx1 = ell + 1;
        m1_total = 0.0;
        m2_total = 0.0;
        set_A_size = 0.0;
        set_B_size = 0.0;
        for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is ik for r.ik; or ij for mu.ij
            if(!R_IsNA(param_vec[loop_idx])) {
                if(c_vec(loop_idx)==ell_idx1) {
                    if(lambda_vec(loop_idx) == 1) { // Set 'A' in the full conditionals
                        m1_total += param_vec(loop_idx);
                        set_A_size++;
                    }
                    else if(lambda_vec[loop_idx] == 0){ // lambda_vec[loop_idx]==0; set 'B' in full conditionals
                        m2_total += ( upsilon - (1-w_vec(ell)) * param_vec(loop_idx) )/w_vec(ell);
                        set_B_size++;
                    }
                    else {
                        Rcpp::Rcout << "loop_idx=" << loop_idx << std::endl;
                        Rcpp::Rcout << "lambda=" << lambda_vec[loop_idx];
                        stop("Inconsistent lambda value in fnSampEtaVec");
                    }
                }
            }
        }
        if(set_A_size > 0) {
            m1 = m1_total/set_A_size;
        }
        else {
            m1 = 0;
        }
        if(set_B_size > 0) {
            m2 = m2_total/set_B_size;
        }
        else {
            m2 = 0;
        }
        rho1 = (1.0/u_2) * set_A_size;
        rho2 = (1.0/u_2) * (w_vec(ell)/(1-w_vec(ell))) * set_B_size;
        norm_precision = rho1 + rho2 + rho3;
        norm_variance = 1.0/norm_precision;
        norm_mean = (rho1*m1 + rho2*m2 + rho3*m3)/norm_precision;

        // Rcpp::Rcout << "m1=" << m1 << std::endl;
        // Rcpp::Rcout << "m2=" << m2 << std::endl;
        // Rcpp::Rcout << "m3=" << m3 << std::endl;
        // Rcpp::Rcout << "norm_mean=" << norm_mean << std::endl;
        // Rcpp::Rcout << "m1_total=" << m1_total << std::endl;
        // Rcpp::Rcout << "m2_total=" << m2_total << std::endl;
        // Rcpp::Rcout << "|A|=" << set_A_size << std::endl;
        // Rcpp::Rcout << "|B|=" << set_B_size << std::endl;
        // Rcpp::Rcout << "rho1=" << rho1 << std::endl;
        // Rcpp::Rcout << "rho2=" << rho2 << std::endl;
        // Rcpp::Rcout << "rho3=" << rho3 << std::endl;

        eta_vec(ell) = R::rnorm(norm_mean, sqrt(norm_variance));
    }

    return eta_vec;
}

double log_sum_exp(double u, double v) {
    // See http://andrewgelman.com/2016/06/11/log-sum-of-exponentials/
    return std::max(u, v) + log(exp(u - std::max(u, v)) + exp(v - std::max(u, v)));
}

// [[Rcpp::export]]
NumericVector fnLogMeanConstraintPrior(const double & param, NumericVector psi_vec, NumericVector eta_vec, NumericVector w_vec, const double & upsilon, const double & u_2) {
    // Lambda has been integrated out in this formulation
    const double & L = psi_vec.size();
    NumericVector log_probs(L);

    double phi1;
    double phi2;
    double mean1;
    double mean2;
    for(int ell=0; ell<L; ell++) {
        mean1 = eta_vec(ell);
        mean2 = (upsilon - w_vec(ell)*eta_vec(ell))/(1-w_vec(ell));
        phi1 = R::dnorm(param, mean1, sqrt(u_2), TRUE);
        phi2 = R::dnorm(param, mean2, sqrt(u_2), TRUE);
        log_probs(ell) = log(psi_vec(ell)) + log_sum_exp(log(w_vec(ell)) + phi1, log(1-w_vec(ell)) + phi2);
    }

    return log_probs;
}

// [[Rcpp::export]]
int fnSampC(const double & param, NumericVector psi_vec,
            NumericVector w_vec, NumericVector eta_vec,
            const double & upsilon, const double & u_2) {
    const double & L = psi_vec.size();

    NumericVector log_probs = fnLogMeanConstraintPrior(param, psi_vec, eta_vec, w_vec, upsilon, u_2);
    NumericVector probs = exp(log_probs - max(log_probs));

    return sample(L, 1, FALSE, probs)(0);
}

// [[Rcpp::export]]
IntegerVector fnSampCVec(IntegerVector c_vec, IntegerVector d_vec,
                NumericVector param_vec, NumericVector psi_vec,
                NumericVector w_vec, NumericVector eta_vec,
                const double & upsilon, const double & u_2) {

    int new_c;
    int old_c;
    for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is ik for r; or ij for w_\ell^\theta
        if(!R_IsNA(param_vec[loop_idx])) {
            old_c = c_vec(loop_idx);
            new_c = fnSampC(param_vec(loop_idx), psi_vec, w_vec, eta_vec, upsilon, u_2);

            c_vec(loop_idx) = new_c;

            if(new_c != old_c) { // Need to update cluster counts
                d_vec(old_c-1) -= 1; // Don't forget about C's 0 indexing here
                d_vec(new_c-1) += 1;
            }
        }
    }

    return c_vec;
}

// [[Rcpp::export]]
IntegerVector fnSampLambdaVec(IntegerVector lambda_vec, NumericVector param_vec,
                     IntegerVector c_vec, NumericVector w_vec,
                     NumericVector eta_vec, const double & upsilon,
                     const double & u_2) {
    int new_lambda;
    NumericVector log_probs(2);
    NumericVector probs(2);
    double ell_idx0;
    double ell_idx1;
    double w_ell;
    double eta_ell;
    double mean1;
    double mean2;
    double p; // Probability that lambda=1
    for(int loop_idx=0; loop_idx < c_vec.size(); loop_idx++) { // Loop index is tk for w_\ell^r; or j for w_\ell^\theta
        if(!R_IsNA(param_vec[loop_idx])) {
            ell_idx1 = c_vec(loop_idx);
            ell_idx0 = ell_idx1 - 1;
            w_ell = w_vec(ell_idx0);
            eta_ell = eta_vec(ell_idx0);

            mean1 = eta_ell;
            mean2 = (upsilon - w_ell*eta_ell)/(1-w_ell);

            log_probs(1) = log(w_ell) + R::dnorm(param_vec(loop_idx), mean1, sqrt(u_2), TRUE);
            log_probs(0) = log(1-w_ell) + R::dnorm(param_vec(loop_idx), mean2, sqrt(u_2), TRUE);

            probs = exp(log_probs - max(log_probs));
            p = probs(1)/sum(probs);
            new_lambda = R::rbinom(1, p);

            lambda_vec(loop_idx) = new_lambda;

            // Rcpp::Rcout << "log_probs=" << log_probs << std::endl;
            // Rcpp::Rcout << "probs=" << probs << std::endl;
            // Rcpp::Rcout << "p=" << p << std::endl;
        }
    }

    return lambda_vec;
}

#endif
