#ifndef UTILITIES_H
#define UTILITIES_H

#include <Rcpp.h>
using namespace Rcpp;

/* Functions ################################################################ */
// [[Rcpp::export]]
bool fnAprxEql(double x, double y) {
    const double epsilon = 1.5E-7;
    const double tolerance = 1.5E-7;

    bool eql = false;
    if((fabs(x) < tolerance) && (fabs(y) < tolerance)) {
        eql = true;
    }
    else {
        eql = fabs(x-y)/(fabs(x) + fabs(y) + epsilon) < tolerance;
    }

    return(eql);
}

// [[Rcpp::export]]
bool fnAcceptProposal(const double & curr_log_lik, const double & prop_log_lik, const double log_curr_to_prop_prob = 0, const double log_prop_to_curr_prob = 0) {
    double u = runif(1)(0);
    bool accept = FALSE;
    if(log(u) <= (prop_log_lik - curr_log_lik + log_prop_to_curr_prob - log_curr_to_prop_prob)) {
        accept = TRUE;
    }

    return accept;
}

// [[Rcpp::export]]
void fnMuMatCorrect(NumericMatrix Mu_mat, NumericMatrix Delta_mat, NumericVector r_vec, NumericMatrix Alpha_mat, NumericMatrix Theta_mat, NumericVector rep_K, NumericVector subject, NumericVector condition, Nullable<NumericMatrix> Xtbeta_mat_=R_NilValue) {
    // const int & n = Alpha_mat.nrow();
    const int & J = Alpha_mat.ncol();
    const int & N= Mu_mat.nrow();

    NumericMatrix Xtbeta_mat(1, 1);
    double log_calc_mu_ikj;
    double calc_mu_ikj;
    int k = -1;
    int i = -1;
    for(int j=0; j<J; j++) {
        for(int ik=0; ik<N; ik++) {
            k = condition[ik] - 1;
            i = subject[ik] - 1;
            // if( (Delta_mat(ik, j) == 1) && !fnAprxEql(log(Mu_mat(ik, j)), r_vec[ik] + Alpha_mat(i, j) + Theta_mat(k, j)) ) {
            if(Xtbeta_mat_.isNull()) {
                log_calc_mu_ikj = r_vec[ik] + Alpha_mat(i, j) + Theta_mat(k, j);
            }
            else {
                Xtbeta_mat = NumericMatrix(Xtbeta_mat_);
                log_calc_mu_ikj = r_vec[ik] + Alpha_mat(i, j) + Theta_mat(k, j) + Xtbeta_mat(ik, j);
            }
            calc_mu_ikj = exp(log_calc_mu_ikj);
            if( (Delta_mat(ik, j) == 1) && !fnAprxEql(Mu_mat(ik, j), calc_mu_ikj) ){
                Rcpp::Rcout << "i=" << i << std::endl;
                Rcpp::Rcout << "k=" << k << std::endl;
                Rcpp::Rcout << "ik=" << ik << std::endl;
                Rcpp::Rcout << "j=" << j << std::endl;
                Rcpp::Rcout << "r=" << r_vec[ik] << std::endl;
                Rcpp::Rcout << "alpha=" << Alpha_mat(i, j) << std::endl;
                Rcpp::Rcout << "theta=" << Theta_mat(k, j) << std::endl;
                Rcpp::Rcout <<  std::setprecision(30) << "Mu.ikj=" << Mu_mat(ik, j) << std::endl;
                Rcpp::Rcout <<  std::setprecision(30) << "Calculated Mu.ikj=" << calc_mu_ikj << std::endl;
                Rcpp::Rcout <<  std::setprecision(30) << "log(Mu.ikj)=" << log(Mu_mat(ik, j)) << std::endl;
                if(Xtbeta_mat_.isNotNull()) {
                    Rcpp::Rcout <<  "xtbeta=" << Xtbeta_mat(ik, j) << std::endl;
                }
                Rcpp::Rcout <<  std::setprecision(30) << "Calculated log(Mu.ikj)=" << log_calc_mu_ikj << std::endl << std::endl;
                stop("Mu.ikj inconsistent");
            }
            if( (Delta_mat(ik, j) == 0) && !(Mu_mat(ik, j)==0) ) {
                Rcpp::Rcout << "i=" << i << std::endl;
                Rcpp::Rcout << "k=" << k << std::endl;
                Rcpp::Rcout << "ik=" << ik << std::endl;
                Rcpp::Rcout << "j=" << j << std::endl;
                Rcpp::Rcout << "r=" << r_vec[ik] << std::endl;
                Rcpp::Rcout << "alpha=" << Alpha_mat(i, j) << std::endl;
                Rcpp::Rcout << "theta=" << Theta_mat(k, j) << std::endl;
                Rcpp::Rcout << "Mu.ikj=" << std::setprecision(30) << Mu_mat(ik, j) << std::endl;
                stop("delta.ikj = 0 but mu.ikj != 0");
            }
        }
    }
}

#endif
