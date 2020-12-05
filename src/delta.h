#ifndef DELTA_H
#define DELTA_H

#include <Rcpp.h>
#include "likelihood.h"
#include "utilities.h"
using namespace Rcpp;

struct alpha_draw_strct {
    int c_ij;
    int lambda_ij;
    double alpha_ij;
};

struct theta_draw_strct {
    double c_kj;
    double theta_kj;
};

alpha_draw_strct fnDrawAlphaij(NumericVector psi_alpha_vec, NumericVector w_alpha_vec, NumericVector eta_alpha_vec, const double & upsilon_alpha, const double & u_alpha_2, const int & a_ij) {
    alpha_draw_strct alpha_draw;
    if(a_ij == 0) {
        alpha_draw.c_ij = NA_INTEGER;
        alpha_draw.lambda_ij = NA_INTEGER;
        alpha_draw.alpha_ij = NA_REAL;
    }
    else if(a_ij == 1) {
        const int & L = psi_alpha_vec.size();

        alpha_draw.c_ij = sample(L, 1, FALSE, psi_alpha_vec)(0);
        int l_idx0 = alpha_draw.c_ij - 1;
        alpha_draw.lambda_ij = R::rbinom(1.0, w_alpha_vec[l_idx0]);

        if(alpha_draw.lambda_ij==0) {
            alpha_draw.alpha_ij = R::rnorm((upsilon_alpha - w_alpha_vec[l_idx0] * eta_alpha_vec[l_idx0])/(1.0 - w_alpha_vec[l_idx0]), sqrt(u_alpha_2));
        }
        else if(alpha_draw.lambda_ij==1) {
            alpha_draw.alpha_ij = R::rnorm(eta_alpha_vec[l_idx0], sqrt(u_alpha_2));
        }
        else {
            stop("Unexpected lambda value in fnDrawAlphaij");
        }
    }
    else {
        stop("Unexpected a_ij value passed to fnDrawAlphaij");
    }


    return alpha_draw;
}

theta_draw_strct fnDrawThetakj(NumericVector psi_theta_vec, NumericVector theta_star_k_vec, const double & sigma_2_theta_k, const int & b_kj, const int & delta_kj_sum) {
    theta_draw_strct theta_draw;
    if(b_kj==0) {
        theta_draw.c_kj = NA_REAL;
        theta_draw.theta_kj = NA_REAL;
    }
    else if(b_kj==1 && delta_kj_sum==0) {
        theta_draw.c_kj = NA_REAL;
        theta_draw.theta_kj = 0.0;
    }
    else {
        const int & L = psi_theta_vec.size();

        theta_draw.c_kj = sample(L, 1, FALSE, psi_theta_vec)(0);
        int l_idx0 = theta_draw.c_kj - 1;
        theta_draw.theta_kj = R::rnorm(theta_star_k_vec[l_idx0], sqrt(sigma_2_theta_k));
    }

    return theta_draw;
}

int fnGetAij(const int & i_idx0, const int & j_idx0, NumericMatrix Delta_mat, NumericVector subject) {
    // Could be improved by taking advantage of the fact that the subjects are in order
    const int & N = Delta_mat.nrow();

    int a_ij=0;
    for(int ik=0; ik<N; ik++) {
        if( ((subject[ik]-1)==i_idx0) && (Delta_mat(ik, j_idx0)==1) ) {
            a_ij=1;
            break;
        }
    }

    return a_ij;
}

int fnGetBkj(const int & k_idx0, const int & j_idx0, NumericMatrix Delta_mat, NumericVector condition) {
    const int & N = Delta_mat.nrow();

    int b_kj=0;
    for(int ik=0; ik<N; ik++) {
        if( ((condition[ik]-1)==k_idx0) && Delta_mat(ik, j_idx0)==1) {
            b_kj=1;
            break;
        }
    }

    return b_kj;
}

int fnDeltakjSum(const int & k_idx0, const int & j_idx0, NumericMatrix Delta_mat, NumericVector condition, const double & K, int *k_prime) {
    int b_kj;
    int delta_kj_sum=0;
    for(int k_iter=0; k_iter<K; k_iter++) {
        if(k_iter!=k_idx0) {
            b_kj = fnGetBkj(k_iter, j_idx0, Delta_mat, condition);
            if(b_kj==1) *k_prime = k_iter;
            delta_kj_sum += b_kj;
        }
    }
    if(delta_kj_sum!=1) *k_prime = NA_INTEGER;

    return delta_kj_sum;
}

NumericVector fnAlphajVecRep(const double & i, const double & alpha_ij_prop, NumericVector alpha_j_vec_curr, NumericVector subject, const double & N) {
    int i_idx0;
    NumericVector alpha_j_vec_rep(N);
    for(int ik=0; ik<N; ik++) {
        i_idx0 = subject[ik] - 1;
        if(i==i_idx0) {
            alpha_j_vec_rep[ik] = alpha_ij_prop;
        }
        else { // Some other subject
            alpha_j_vec_rep[ik] = alpha_j_vec_curr[i_idx0];
        }
    }

    return(alpha_j_vec_rep);
}

NumericVector fnThetajVecRep(const double & k, const double & k_prime, const double & theta_kj_prop, const double & theta_kj_prime_prop, NumericVector theta_j_vec_curr, bool theta_prime_change, NumericVector condition, const double & N) {
    double k_idx0;
    NumericVector theta_j_vec_rep(N);
    for(int ik=0; ik<N; ik++) {
        k_idx0 = condition[ik] - 1;
        if(k == k_idx0) {
            theta_j_vec_rep[ik] = theta_kj_prop;
        }
        else if(theta_prime_change && (k_idx0==k_prime)) {
            theta_j_vec_rep[ik] = theta_kj_prime_prop;
        }
        else { // Some other condition
            theta_j_vec_rep[ik] = theta_j_vec_curr[k_idx0];
        }
    }

    return(theta_j_vec_rep);
}

void fnSampDelta(const int & ik, const int & i_idx0, const int & j_idx0,
                 const int & k_idx0,  const int & delta_kj_sum,
                 const int & k_prime, NumericVector psi_alpha_vec,
                 NumericVector w_alpha_vec, NumericVector eta_alpha_vec, IntegerVector c_alpha_vec, IntegerVector d_alpha_vec, IntegerVector lambda_alpha_vec, NumericMatrix Alpha_mat,
                 NumericVector psi_theta_vec, NumericMatrix Theta_star_mat, NumericMatrix Theta_mat, NumericMatrix C_theta_mat, IntegerVector d_theta_vec,
                 NumericVector sigma_2_theta_vec, NumericMatrix Delta_mat, NumericVector r_vec, NumericVector s_vec, NumericMatrix Mu_mat, NumericMatrix Epsilon_mat, NumericMatrix Y_mat,
                 NumericVector subject, NumericVector condition, NumericVector rep_K, NumericMatrix Xtbeta_mat,
                 const double & upsilon_alpha, const double & u_alpha_2) {
    const double & N = Mu_mat.nrow();
    const int & n = Alpha_mat.nrow();
    const int ij = j_idx0 * n + i_idx0;

    const double delta_ikj_curr = Delta_mat(ik, j_idx0);
    const double delta_ikj_prop = 1.0 - delta_ikj_curr;

    // Determine a_ij_0 and b_kj_0
    int a_ij_0 = fnGetAij(i_idx0, j_idx0, Delta_mat, subject);
    int b_kj_0 = fnGetBkj(k_idx0, j_idx0, Delta_mat, condition);

    // Determine a_ij_1 and b_kj_1
    Delta_mat(ik, j_idx0) = delta_ikj_prop; // Set temporarily to get a.ij and b.kj
    int a_ij_1 = fnGetAij(i_idx0, j_idx0, Delta_mat, subject);
    int b_kj_1 = fnGetBkj(k_idx0, j_idx0, Delta_mat, condition);
    Delta_mat(ik, j_idx0) = delta_ikj_curr;

    bool alpha_change = a_ij_0 != a_ij_1;
    bool theta_change = b_kj_0 != b_kj_1;
    bool theta_prime_change = theta_change && (delta_kj_sum==1);
    if(theta_prime_change && IntegerVector::is_na(k_prime)) stop("k_prime is NA when theta_prime_change=1");

    int lambda_ij_prop;
    int c_ij_prop;
    double alpha_ij_prop;
    double c_kj_prop;
    double theta_kj_prop;
    double c_kj_prime_prop;
    double theta_kj_prime_prop;
    // Draw alpha_ij if necessary
    if(alpha_change) {
        alpha_draw_strct alpha_draw = fnDrawAlphaij(psi_alpha_vec, w_alpha_vec, eta_alpha_vec, upsilon_alpha, u_alpha_2, a_ij_1);
        alpha_ij_prop = alpha_draw.alpha_ij;
        c_ij_prop = alpha_draw.c_ij;
        lambda_ij_prop = alpha_draw.lambda_ij;
    }
    else {
        alpha_ij_prop = Alpha_mat(i_idx0, j_idx0);
        c_ij_prop = c_alpha_vec[ij];
        lambda_ij_prop = lambda_alpha_vec[ij];
    }

    // Draw theta_kj if necessary
    if(theta_change) {
        theta_draw_strct theta_draw = fnDrawThetakj(psi_theta_vec, Theta_star_mat(k_idx0, _), sigma_2_theta_vec[k_idx0], b_kj_1, delta_kj_sum);
        theta_kj_prop = theta_draw.theta_kj;
        c_kj_prop = theta_draw.c_kj;
    }
    else {
        theta_kj_prop = Theta_mat(k_idx0, j_idx0);
        c_kj_prop = C_theta_mat(k_idx0, j_idx0);
    }

    // Draw theta_kprime_j if necessary
    if(theta_prime_change) {
        const int b_kj_prime = 1.0;
        double delta_kj_sum_prime;
        if(delta_ikj_curr == 0) {
            delta_kj_sum_prime = 1.0;
        }
        else { // delta_ikj_curr == 1
            delta_kj_sum_prime = 0.0;
        }
        theta_draw_strct theta_prime_draw = fnDrawThetakj(psi_theta_vec, Theta_star_mat(k_prime, _), sigma_2_theta_vec[k_prime], b_kj_prime, delta_kj_sum_prime);
        theta_kj_prime_prop = theta_prime_draw.theta_kj;
        c_kj_prime_prop = theta_prime_draw.c_kj;
    }
    else {
        // Set to some nonsense values, can be removed
        theta_kj_prime_prop = -100;
        c_kj_prime_prop = -100;
    }

    NumericVector alpha_j_vec_rep = fnAlphajVecRep(i_idx0, alpha_ij_prop, Alpha_mat(_, j_idx0), subject, N);
    NumericVector theta_j_vec_rep = fnThetajVecRep(k_idx0, k_prime, theta_kj_prop, theta_kj_prime_prop, Theta_mat(_, j_idx0), theta_prime_change, condition, N);
    NumericVector delta_j_vec_prop(N);
    for(int ik_iter=0; ik_iter<N; ik_iter++) delta_j_vec_prop[ik_iter] = Delta_mat(ik_iter, j_idx0);
    delta_j_vec_prop[ik] = delta_ikj_prop;
    NumericVector mu_j_vec_prop = exp(r_vec + alpha_j_vec_rep + theta_j_vec_rep + Xtbeta_mat(_, j_idx0));
    // for(int ik_iter=0; ik_iter<N; ik_iter++) if(R_IsNA(mu_j_vec_prop[ik_iter]) || (delta_j_vec_prop[ik_iter]==0)) mu_j_vec_prop[ik_iter] = 0.0;
    for(int ik_iter=0; ik_iter<N; ik_iter++) if(delta_j_vec_prop[ik_iter]==0) mu_j_vec_prop[ik_iter] = 0.0;

    double prior_part_curr = (delta_ikj_curr == 1) ? log1p(-Epsilon_mat(k_idx0, j_idx0)) : log(Epsilon_mat(k_idx0, j_idx0));
    double prior_part_prop = (delta_ikj_prop == 1) ? log1p(-Epsilon_mat(k_idx0, j_idx0)) : log(Epsilon_mat(k_idx0, j_idx0));
    double lik_part_curr = fnOtuLogLik(s_vec[j_idx0], Mu_mat(_, j_idx0), rep_K, Y_mat(_, j_idx0), Delta_mat(_, j_idx0));
    double lik_part_prop = fnOtuLogLik(s_vec[j_idx0], mu_j_vec_prop, rep_K, Y_mat(_, j_idx0), delta_j_vec_prop);
    double curr_log_lik = prior_part_curr + lik_part_curr;
    double prop_log_lik = prior_part_prop + lik_part_prop;

    if(fnAcceptProposal(curr_log_lik, prop_log_lik)) {
        if(!IntegerVector::is_na(c_alpha_vec[ij])) d_alpha_vec[c_alpha_vec[ij]-1] -= 1;
        c_alpha_vec[ij] = c_ij_prop;
        if(!IntegerVector::is_na(c_alpha_vec[ij])) d_alpha_vec[c_alpha_vec[ij]-1] += 1;
        lambda_alpha_vec[ij] = lambda_ij_prop;
        Alpha_mat(i_idx0, j_idx0) = alpha_ij_prop;

        if(!IntegerVector::is_na(C_theta_mat(k_idx0, j_idx0))) d_theta_vec[C_theta_mat(k_idx0, j_idx0) - 1] -= 1;
        C_theta_mat(k_idx0, j_idx0) = c_kj_prop;
        if(!IntegerVector::is_na(C_theta_mat(k_idx0, j_idx0))) d_theta_vec[C_theta_mat(k_idx0, j_idx0) - 1] += 1;
        Theta_mat(k_idx0, j_idx0) = theta_kj_prop;
        if(theta_prime_change) {
            if(!IntegerVector::is_na(C_theta_mat(k_prime, j_idx0))) d_theta_vec[C_theta_mat(k_prime, j_idx0) - 1] -= 1;
            C_theta_mat(k_prime, j_idx0) = c_kj_prime_prop;
            if(!IntegerVector::is_na(C_theta_mat(k_prime, j_idx0))) d_theta_vec[C_theta_mat(k_prime, j_idx0) - 1] += 1;
            Theta_mat(k_prime, j_idx0) = theta_kj_prime_prop;
        }

        Mu_mat(_, j_idx0) = mu_j_vec_prop;
        Delta_mat(ik, j_idx0) = delta_ikj_prop;
    }
}

// [[Rcpp::export]]
NumericMatrix fnSampDeltaMat(NumericMatrix Delta_mat, NumericMatrix Epsilon_mat, NumericMatrix Theta_mat, NumericMatrix C_theta_mat, IntegerVector d_theta_vec, NumericMatrix Theta_star_mat, NumericVector psi_theta_vec, NumericVector sigma_2_theta_vec, NumericMatrix Alpha_mat, NumericVector r_vec, NumericVector s_vec, NumericMatrix Mu_mat, NumericMatrix Y_mat, NumericVector condition, NumericVector subject, const double & theta_kj_proposal_sd,
                             NumericVector psi_alpha_vec, NumericVector w_alpha_vec, NumericVector eta_alpha_vec, IntegerVector c_alpha_vec, IntegerVector d_alpha_vec, IntegerVector lambda_alpha_vec, NumericVector rep_K, NumericMatrix Xtbeta_mat, const double & upsilon_alpha, const double & u_alpha_2) {
    const int & N = Delta_mat.nrow();
    // const int & n = Alpha_mat.nrow();
    const int & J = Delta_mat.ncol();
    const int & K = Theta_mat.nrow();

    NumericVector mu_j_vec_prop(N);
    int i_idx0;
    int k_idx0;
    int delta_kj_sum;
    int k_prime;
    for(int j=0; j<J; j++) {
        for(int ik=0; ik<N; ik++) {
            if(Y_mat(ik, j) == 0) { // If Y.ikj>0, delta.ikj=1
                i_idx0 = subject[ik] - 1;
                k_idx0 = condition[ik] - 1;
                delta_kj_sum = fnDeltakjSum(k_idx0, j, Delta_mat, condition, K, &k_prime);
                fnSampDelta(ik, i_idx0, j, k_idx0, delta_kj_sum,
                            k_prime, psi_alpha_vec, w_alpha_vec,
                            eta_alpha_vec, c_alpha_vec, d_alpha_vec, lambda_alpha_vec, Alpha_mat, psi_theta_vec, Theta_star_mat, Theta_mat, C_theta_mat, d_theta_vec,
                            sigma_2_theta_vec,  Delta_mat, r_vec, s_vec, Mu_mat, Epsilon_mat, Y_mat, subject,
                            condition, rep_K, Xtbeta_mat, upsilon_alpha,  u_alpha_2);
            }
        }
    }

    return(Delta_mat);
}

#endif
