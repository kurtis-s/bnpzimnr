#' Load everything in INIT into the parent frame
#'
#' @export
fnInitParams <- function(INIT) {
    for(param.name in names(INIT)) {
        assign(param.name, INIT[[param.name]], envir=parent.frame())
    }
}

#' Fit BNP-ZIMNR
#'
#' Used to fit Bayesian nonparametric multivariate NB regression model with
#' zero-inflation (BNP-ZIMNR). See the accompanying paper in the references for
#' a more detailed description of the model and its implementation.
#'
#' @param dat data frame of count data with samples in rows and counts in
#'   columns.  The first two columns of the dataframe should have column names
#'   "subject" and "condition", indicating the grouping factor and sample
#'   condition under which the sample was taken. The model implementation
#'   expects "subject" and "condition" to be ordered numeric labels, starting
#'   with an index of 1 for the first subject.  See the dataset included in the
#'   package for an example.  The remaining columns are compromised of
#'   non-negative integer counts, such as OTU counts in a microbiome study.
#' @param n.save number of samples to take.  If \code{thin}=1, the number of
#'   iterations to run the MCMC chain
#' @param thin thinning interval
#' @param INIT list of parameters used to initialize the MCMC chain with elements:
#' \itemize{
#'     \item s.vec; Overdispersion parameters
#'     \item r.vec; Normalization parameters
#'     \item c.r.vec; Cluster assignments for normalization parameters, 1-indexed
#'     \item eta.r.vec; Mean parameters of the first component of the mixture-of-mixtures prior on the normalization parameters
#'     \item w.r.vec; Mixture weights for the two Gaussian components of the normalization parameters
#'     \item lambda.r.vec; Indicators for the Gaussian component of the normalization parameters
#'     \item Delta.mat; Zero inflation indicators
#'     \item Epsilon.mat; Zero inflation probabilities (Should correspond with Xi.mat)
#'     \item Xi.mat; Probit of Epsilon.mat
#'     \item C.xi.mat; Cluster assignments for zero inflation probabilities, 1-indexed
#'     \item Xi.star.mat; Cluster means for zero inflation probabilities
#'     \item sigma.2.xi.vec; Variances of Xi.mat mixture distribution
#'     \item Alpha.mat; Baseline abundances
#'     \item c.alpha.vec; Cluster assignments for baseline abundances
#'     \item eta.alpha.vec; Mean parameters of the first component of the mixture-of-mixtures prior on the baseline abundances
#'     \item w.alpha.vec;  Mixture weights for the two Gaussian components of the baseline abundances
#'     \item lambda.alpha.vec; Indicators for the Gaussian component of the baseline abundances
#'     \item Theta.mat; Change in abundance compared to the baseline abundance
#'     \item C.theta.mat; Cluster assignments for Theta.mat
#'     \item Theta.star.mat; Mean parameters for Theta.mat mixture
#'     \item sigma.2.theta.vec; Variance parameters for Theta.mat mixture
#' }
#' @param PRIOR list of hyperparameters with elements:
#' \itemize{
#'     \item a.s; mean of s.j prior
#'     \item b.s.2; variance of s.j prior
#'     \item a.psi.r; dirichlet concentration parameter for psi.r
#'     \item a.w.r; theta shape parameter for w.r
#'     \item b.w.r; theta shape parameter for w.r
#'     \item u.r.2; normal variance for r.ik mixture
#'     \item upsilon.r; r.ik mean constraint
#'     \item b.eta.r.2; normal variance for eta.r
#'     \item a.psi.alpha; Dirichlet concentration parameter for psi.alpha
#'     \item a.w.alpha; theta shape parameter for w.alpha
#'     \item b.w.alpha; theta shape parameter for w.alpha
#'     \item u.alpha.2; normal variance for alpha mixture
#'     \item upsilon.alpha; alpha mean constraint
#'     \item b.eta.alpha.2; normal variance for eta.alpha
#'     \item rho.theta; Dirichlet-Process concentration parameter for theta
#'     \item a.theta.star; normal mean for theta.star
#'     \item b.theta.star.2; normal variance for theta.star
#'     \item a.sigma.theta; shape parameter for theta variance
#'     \item b.sigma.theta; scale parameter for theta variance
#'     \item rho.xi; Dirichlet-Process concentration parameter for xi
#'     \item a.xi.star; mean for xi.star
#'     \item b.xi.star.2; variance for xi.star
#'     \item a.sigma.xi; shape parameter for xi variance
#'     \item b.sigma.xi; scale parameter for xi variance
#' }
#' @param SETTINGS list of MCMC run settings with elements:
#' \itemize{
#'     \item print.mcmc; Print MCMC progress after print.mcmc iterations
#'     \item s.j.tilde.proposal.sd; Standard deviation for s.j Metropolis proposal
#'     \item L.r.trunc; Truncation level for r mixture
#'     \item w.ell.r.proposal.sd; Standard deviation for w.ell.r Metropolis proposal
#'     \item r.ik.proposal.sd; Standard deviation for r Metropolis proposal
#'     \item L.alpha.trunc; Truncation level for alpha mixture
#'     \item w.ell.alpha.proposal.sd; Standard deviation for w.ell.alpha Metropolis proposal
#'     \item alpha.ij.proposal.sd; Standard deviation for alpha Metropolis proposal
#'     \item L.theta.trunc; Truncation level for theta mixture
#'     \item theta.kj.proposal.sd; Standard deviation for theta Metropolis proposal
#'     \item L.xi.trunc; Truncation level for xi mixture
#'     \item xi.kj.proposal.sd; Standard deviation for xi Metropolis proposal
#' }
#' @param adapt.stop for the adaptive Metropolis samplers stop adapting
#'   Metropolis proposal SDs after adapt.stop iterations (for example, stop
#'   adapting after some burn-in period).  If set to NA the proposal SDs will
#'   continued to be modified throughout the entire MCMC.
#' @param sanity.checks boolean for whether to periodically check of the
#'     parameters are in a consistent state.
#' @return A list whose elements contain samples of the model paramters; one
#'   sample per row
#' @references \emph{A Bayesian Nonparametric Analysis for Zero Inflated
#'   Multivariate Count Data with Application to Microbiome Study}
#' @example examples/example-zimnr.R
#' @export
zimnr <- function(dat, n.save, thin, INIT, PRIOR, SETTINGS, adapt.stop=NA, delta.fixed=FALSE, poisson.model=FALSE, fixed.norm=FALSE, X.mat=NA) {
    INIT$d.r.vec <- as.integer(INIT$d.r.vec) # Modified in place, make sure it's the correct type
    INIT$d.alpha.vec <- as.integer(INIT$d.alpha.vec) # Modified in place, make sure it's the correct type
    d.theta.vec <- as.integer(tabulate(INIT$C.theta.mat, nbins=SETTINGS$L.theta.trunc))
    d.theta.cumsum.vec <- as.integer(rev(cumsum(rev(d.theta.vec)))[-1])
    d.xi.vec <- as.integer(tabulate(INIT$C.xi.mat, nbins=SETTINGS$L.xi.trunc))
    d.xi.cumsum.vec <- as.integer(rev(cumsum(rev(d.xi.vec)))[-1])
    INIT$Epsilon.mat <- matrix(as.numeric(INIT$Epsilon.mat), nrow=nrow(INIT$Epsilon.mat), ncol=ncol(INIT$Epsilon.mat)) # Modified in place, make sure it's the correct type
    # INIT$Mu.mat <- matrix(as.numeric(INIT$Mu.mat), nrow=nrow(INIT$Mu.mat), ncol=ncol(INIT$Mu.mat)) # Modified in place, make sure it's the correct type
    fnInitParams(INIT)
    Y.mat <- as.matrix(dat[,-c(1:2)])
    # Y.kj.prod.mat <- apply(Y.mat, 2, function(otu) tapply(otu==0, factor(dat$condition), FUN=prod))

    if(poisson.model) s.vec <- rep(0, ncol(Y.mat))
    if(fixed.norm) {
        r.vec <- log(rowSums(Y.mat))
    }

    # Make sure Mu.mat is set to the correct values initially
    if(length(X.mat)==1 && is.na(X.mat)) {
        Beta.mat <- NA
    }
    Mu.mat <- fnGetMuMat(r.vec, Alpha.mat, Theta.mat, Delta.mat, s.vec, dat$subject, dat$condition, X.mat, Beta.mat)

    # The delta sampler expects 'condition' and 'subject' to be ordered and in a particular format
    # if(!isTRUE(all.equal(dat$subject, rep(1:length(unique(dat$subject)), each=length(unique(dat$condition)))))) stop("subject in unexpected order or has unexpected labels")
    # Check of the subjects are in order
    if(!isTRUE(all(dat$subject == sort(dat$subject)))) stop("expected subject observations to be grouped together in order")
    if(!isTRUE(all.equal(unique(dat$subject), 1:length(unique(dat$subject))))) stop("subject in unexpected order or has unexpected labels")
    # if(!isTRUE(all.equal(dat$condition, rep(1:length(unique(dat$condition)), times=length(unique(dat$subject)))))) stop("condition in unexpected order or has unexpected labels")
    # Check of conditions are in order for each subject
    if(!isTRUE(all(by(dat$condition, dat$subject, FUN=function(x) all(x==sort(x)))))) stop("condition in unexpected order")

    rep.K <- as.vector(table(dat$subject))

    fnSampRVecAdapt <- overture::Amwg(fnSampRVec, rep(SETTINGS$r.ik.proposal.sd, length(r.vec)), stop.after=adapt.stop)
    fnSampSVecAdapt <- overture::Amwg(fnSampSVec, rep(SETTINGS$s.j.tilde.proposal.sd, length(s.vec)), stop.after=adapt.stop)
    fnSampWRVecAdapt <- overture::Amwg(fnSampWVec, rep(SETTINGS$w.ell.r.proposal.sd, SETTINGS$L.r.trunc), stop.after=adapt.stop)
    # fnSampBetaMatAdapt <- overture::Amwg(fnSampBetaMat, SETTINGS$beta.vec.proposal.sd, stop.after=adapt.stop)
    # fnSampAlphaMatAdapt <- overture::Amwg(fnSampAlphaMat, rep(SETTINGS$alpha.ij.proposal.sd, length(Alpha.mat)), stop.after=adapt.stop)
    # fnSampWAlphaMatAdapt <- overture::Amwg(fnSampWVec, rep(SETTINGS$w.ell.alpha.proposal.sd, SETTINGS$L.alpha.trunc), stop.after=adapt.stop)

    exclude <- NULL

    z <- 0
    fnMcmc <- overture::InitMcmc(n.save=n.save, thin=thin, exclude=exclude)
    samples <- fnMcmc({
        z <- z + 1

        ## s.j
        if(poisson.model) {
            s.vec <- rep(0, ncol(Y.mat))
        } else {
            s.vec <- fnSampSVecAdapt(s.vec, Mu.mat, rep.K, Y.mat, Delta.mat, PRIOR$a.s, PRIOR$b.s.2)
        }

        ## r.ik
        if(fixed.norm) {
            r.vec <- log(rowSums(Y.mat))
        } else {
            psi.r.vec <- fnSampPsi(d.r.vec, PRIOR$a.psi.r)
            w.r.vec <- fnSampWRVecAdapt(w.r.vec, r.vec, lambda.r.vec,
                       c.r.vec, eta.r.vec, PRIOR$upsilon.r, PRIOR$a.w.r,
                       PRIOR$b.w.r, PRIOR$u.r.2)
            eta.r.vec <- fnSampEtaVec(eta.r.vec, r.vec, w.r.vec, lambda.r.vec, c.r.vec,
                         PRIOR$u.r.2, PRIOR$upsilon.r, PRIOR$b.eta.r.2)
            c.r.vec <- fnSampCVec(c.r.vec, d.r.vec, r.vec, psi.r.vec, w.r.vec, eta.r.vec,
                       PRIOR$upsilon.r, PRIOR$u.r.2) # Modifies d.r.vec in place
            lambda.r.vec <- fnSampLambdaVec(lambda.r.vec, r.vec, c.r.vec, w.r.vec,
                            eta.r.vec, PRIOR$upsilon.r, PRIOR$u.r.2)
            r.vec <- fnSampRVecAdapt(r.vec, Mu.mat, s.vec, lambda.r.vec, c.r.vec,
                       eta.r.vec, w.r.vec, rep.K, PRIOR$upsilon.r, PRIOR$u.r.2,
                       Y.mat, Delta.mat) # Modifies Mu.mat in place
        }

        ## alpha.ij
        psi.alpha.vec <- fnSampPsi(d.alpha.vec, PRIOR$a.psi.alpha)
        w.alpha.vec <- fnSampWVec(w.alpha.vec, Alpha.mat, lambda.alpha.vec,
                   c.alpha.vec, eta.alpha.vec, PRIOR$upsilon.alpha, PRIOR$a.w.alpha,
                   PRIOR$b.w.alpha, PRIOR$u.alpha.2, rep(SETTINGS$w.ell.alpha.proposal.sd, SETTINGS$L.alpha.trunc))
        eta.alpha.vec <- fnSampEtaVec(eta.alpha.vec, Alpha.mat, w.alpha.vec, lambda.alpha.vec, c.alpha.vec,
                     PRIOR$u.alpha.2, PRIOR$upsilon.alpha, PRIOR$b.eta.alpha.2)
        c.alpha.vec <- fnSampCVec(c.alpha.vec, d.alpha.vec, Alpha.mat, psi.alpha.vec, w.alpha.vec, eta.alpha.vec,
                   PRIOR$upsilon.alpha, PRIOR$u.alpha.2) # Modifies d.alpha.vec in place
        lambda.alpha.vec <- fnSampLambdaVec(lambda.alpha.vec, Alpha.mat, c.alpha.vec, w.alpha.vec,
                        eta.alpha.vec, PRIOR$upsilon.alpha, PRIOR$u.alpha.2)
        Alpha.mat <- fnSampAlphaMat(Alpha.mat, Mu.mat, s.vec, lambda.alpha.vec, c.alpha.vec,
                                    eta.alpha.vec, w.alpha.vec, rep.K, PRIOR$upsilon.alpha,
                                    PRIOR$u.alpha.2, Y.mat, Delta.mat, SETTINGS$alpha.ij.proposal.sd) # Modifies Mu.mat in place

        ## theta.kj
        v.theta.vec <- fnSampVVec(PRIOR$rho.theta, d.theta.vec, d.theta.cumsum.vec)
        psi.theta.vec <- fnStickBreak(v.theta.vec)
        Theta.star.mat <- fnSampStarMat(Theta.mat, C.theta.mat, sigma.2.theta.vec, PRIOR$a.theta.star, PRIOR$b.theta.star.2, SETTINGS$L.theta.trunc)
        C.theta.mat <- fnSampCMat(C.theta.mat, Theta.mat, d.theta.vec, psi.theta.vec, Theta.star.mat, sigma.2.theta.vec) # Modifies d.theta.vec  vec in place
        d.theta.vec <- d.theta.vec
        d.theta.cumsum.vec <- as.integer(rev(cumsum(rev(d.theta.vec)))[-1])
        B.mat <- fnGetBMat(Delta.mat, dat$condition)
        Delta.sum.mat <- fnGetDeltaSumMat(B.mat)
        Theta.mat <- fnSampThetaMat(Theta.mat, Mu.mat, C.theta.mat, Theta.star.mat, Delta.mat, sigma.2.theta.vec, s.vec, r.vec, Alpha.mat, Y.mat, dat$condition, dat$subject, B.mat, Delta.sum.mat, SETTINGS$theta.kj.proposal.sd) # Modifies Mu.mat in place
        Theta.mat <- fnSampThetaMatMarg(Theta.mat, Mu.mat, psi.theta.vec, Theta.star.mat, Delta.mat, sigma.2.theta.vec, s.vec, r.vec, Alpha.mat, Y.mat, dat$condition, dat$subject, B.mat, Delta.sum.mat, SETTINGS$theta.kj.proposal.sd) # Modifies Mu.mat in place

        ## beta.pj
        if(length(X.mat)==1 && is.na(X.mat)) {
            Xtbeta.mat <- matrix(0, nrow=nrow(Mu.mat), ncol=ncol(Mu.mat)) # The delta sampler expects this quantity to be 0 if X.mat is not supplied
        } else {
            Beta.mat <- fnSampBetaMat(Beta.mat, Mu.mat, Delta.mat, s.vec, tau.2.vec, X.mat, Y.mat, rep.K, rep(SETTINGS$beta.pj.proposal.sd, ncol(Beta.mat)))
            Xtbeta.mat <- X.mat %*% Beta.mat
            tau.2.vec <- fnSampTau2Vec(Beta.mat, PRIOR$a.tau, PRIOR$b.tau)
        }

        ## epsilon.kj
        if(!delta.fixed) {
            v.xi.vec <- fnSampVVec(PRIOR$rho.xi, d.xi.vec, d.xi.cumsum.vec)
            psi.xi.vec <- fnStickBreak(v.xi.vec)
            Xi.star.mat <- fnSampStarMat(Xi.mat, C.xi.mat, sigma.2.xi.vec, PRIOR$a.xi.star, PRIOR$b.xi.star.2, SETTINGS$L.xi.trunc)
            C.xi.mat <- fnSampCMat(C.xi.mat, Xi.mat, d.xi.vec, psi.xi.vec, Xi.star.mat, sigma.2.xi.vec) # Modifies d.xi.vec vec in place
            d.xi.vec <- d.xi.vec
            d.xi.cumsum.vec <- as.integer(rev(cumsum(rev(d.xi.vec)))[-1])
            Xi.mat <- fnSampXiMat(Xi.mat, Epsilon.mat, Delta.mat, Xi.star.mat, C.xi.mat, sigma.2.xi.vec, dat$condition, SETTINGS$xi.kj.proposal.sd) # Modifies Epsilon.mat in place
            Xi.mat <- fnSampXiMatMarg(Xi.mat, Epsilon.mat, Delta.mat, Xi.star.mat, psi.xi.vec, sigma.2.xi.vec, dat$condition, SETTINGS$xi.kj.proposal.sd) # Modifies Epsilon.mat in place
            Epsilon.mat <- Epsilon.mat
            # Xi.star.mat <- fnSampStarMatTrunc(Xi.mat, C.xi.mat, sigma.2.xi.vec, PRIOR$a.xi.star, PRIOR$b.xi.star.2, SETTINGS$L.xi.trunc, PRIOR$xi.star.lbound, PRIOR$xi.star.ubound)
        }

        ## delta.kj
        if(!delta.fixed) {
            Delta.mat <- fnSampDeltaMat(Delta.mat, Epsilon.mat, Theta.mat, C.theta.mat, d.theta.vec, Theta.star.mat, psi.theta.vec, sigma.2.theta.vec, Alpha.mat, r.vec, s.vec, Mu.mat, Y.mat, dat$condition, dat$subject, SETTINGS$theta.kj.proposal.sd,
                   psi.alpha.vec, w.alpha.vec, eta.alpha.vec, c.alpha.vec, d.alpha.vec, lambda.alpha.vec, rep.K, Xtbeta.mat, PRIOR$upsilon.alpha, PRIOR$u.alpha.2) # Modifies Theta.mat, C.theta.mat, Mu.mat, Alpha.mat, lambda.alpha.vec, c.alpha,vec in place
            d.theta.cumsum.vec <- as.integer(rev(cumsum(rev(d.theta.vec)))[-1])
        } else {
            Delta.mat <- Delta.mat
        }

        ## sigma.2
        sigma.2.xi.vec <- fnSampSigma2(Xi.mat, Xi.star.mat, C.xi.mat, PRIOR$a.sigma.xi, PRIOR$b.sigma.xi)
        sigma.2.theta.vec <- fnSampSigma2(Theta.mat, Theta.star.mat, C.theta.mat, PRIOR$a.sigma.theta, PRIOR$b.sigma.theta)
        # sigma.2.xi.vec <- rep(sigma.2.xi, length(sigma.2.xi.vec)) # Only one value for sigma.2 currently, but could be expanded to be sigma.2.k
        # sigma.2.theta.vec <- rep(sigma.2.theta, length(sigma.2.theta.vec)) # Only one value for sigma.2 currently, but could be expanded to be sigma.2.k

        ## Mu.ikj
        Mu.mat <- Mu.mat # So that overture saves Mu.mat
        if(z %% 1000 == 0) Mu.mat <- fnGetMuMat(r.vec, Alpha.mat, Theta.mat, Delta.mat, s.vec, dat$subject, dat$condition, X.mat, Beta.mat) # Reset the Mu.mat calculation periodically

        ## Just for safety.  Can be commented out.
        if(length(X.mat)==1 && is.na(X.mat)) {
            fnMuMatCorrect(Mu.mat, Delta.mat, r.vec, Alpha.mat, Theta.mat, rep.K, dat$subject, dat$condition)
        } else {
            fnMuMatCorrect(Mu.mat, Delta.mat, r.vec, Alpha.mat, Theta.mat, rep.K, dat$subject, dat$condition, Xtbeta.mat)
        }

        # r.ik
        if(!fixed.norm) {
            if(!isTRUE(all.equal(tabulate(c.r.vec, nbins=SETTINGS$L.r.trunc), d.r.vec))) stop("r cluster count in inconsistent state")
            if(!isTRUE(all.equal(sum(psi.r.vec), 1))) stop("Psi.r probabilities don't sum to 1")
        }

        # alpha.ij
        if(!isTRUE(all.equal(sum(psi.alpha.vec), 1))) stop("Psi.alpha probabilities don't sum to 1")
        if(!isTRUE(all.equal(tabulate(c.alpha.vec, nbins=SETTINGS$L.alpha.trunc), d.alpha.vec))) stop("alpha cluster count in inconsistent state")
        A.mat.chk <- fnGetAMat(Delta.mat, dat$subject)
        if(!all(is.na(Alpha.mat[A.mat.chk==0]))) stop("Some alpha.ij are not NA when sum_k delta_ikj = 0")
        if(any(is.na(Alpha.mat[A.mat.chk==1]))) stop("Some alpha.ij are NA when sum_k delta_ikj != 0")
        if(!all(is.na(Alpha.mat[is.na(c.alpha.vec)]))) stop("Some alpha.ij are not NA when c.ij are NA")
        if(!all(is.na(c.alpha.vec[is.na(Alpha.mat)]))) stop("Some c.ij are not NA when alpha.ij are NA")
        if(!all(is.na(lambda.alpha.vec[is.na(Alpha.mat)]))) stop("Some lambda.ij are not NA when alpha.ij are NA")
        if(!all(is.na(Alpha.mat[is.na(lambda.alpha.vec)]))) stop("Some alpha.ij are not NA when lambda.ij are NA")
        if(!all(!is.na(Alpha.mat[!is.na(c.alpha.vec)]))) stop("Some alpha.ij are NA when c.ij is not NA")
        if(!all(!is.na(c.alpha.vec[!is.na(Alpha.mat)]))) stop("Some c.ij are NA when alpha.ij is not NA")
        if(!all(!is.na(Alpha.mat[!is.na(lambda.alpha.vec)]))) stop("Some alpha.ij are NA when lambda.ij is not NA")
        if(!all(!is.na(lambda.alpha.vec[!is.na(Alpha.mat)]))) stop("Some lambda.ij are NA when alpha.ij is not NA")

        # theta.kj
        B.mat.chk <- fnGetBMat(Delta.mat, dat$condition)
        Delta.sum.mat.chk <- fnGetDeltaSumMat(B.mat.chk)
        if(!all(is.na(Theta.mat[B.mat.chk==0]))) stop("Some theta.kj are not NA when sum_i delta_ikj = 0")
        if(any(is.na(Theta.mat[B.mat.chk==1]))) stop("Some theta.kj are NA when sum_i delta_ikj != 0")
        if(!all(Theta.mat[Delta.sum.mat.chk==0] == 0)) stop("Some theta.kj are not 0 when sum_{k^prime != k} b_{kj} = 0")
        if(!all(rev(cumsum(rev(d.theta.vec)))[-1] == d.theta.cumsum.vec)) stop("d.theta.vec and d.theta.cumsum.vec in inconsistent state")
        if(!all(d.theta.vec == tabulate(C.theta.mat, nbins=SETTINGS$L.theta.trunc))) stop("d.theta.vec and C.theta.mat in inconsistent state")
        if(!isTRUE(all.equal(sum(psi.theta.vec), 1))) stop("psi.theta.vec doesn't sum up to 1")
        if(!all(is.na(Theta.mat[is.na(C.theta.mat)]) | Theta.mat[is.na(C.theta.mat)] == 0)) stop("Some theta.kj are not 0 or NA when c.kj is NA")
        if(!all(is.na(C.theta.mat[is.na(Theta.mat)]))) stop("Some c.kj are not NA when theta.kj are NA")

        # epsilon.kj
        if(!delta.fixed) {
            if(!all(rev(cumsum(rev(d.xi.vec)))[-1] == d.xi.cumsum.vec)) stop("d.xi.vec and d.xi.cumsum.vec in inconsistent state")
            if(!all(d.xi.vec == tabulate(C.xi.mat, nbins=SETTINGS$L.xi.trunc))) stop("d.xi.vec and C.xi.mat in inconsistent state")
            if(!isTRUE(all.equal(sum(psi.xi.vec), 1))) stop("psi.xi.vec doesn't sum up to 1")
            if(!isTRUE(all.equal(Epsilon.mat, pnorm(Xi.mat)))) stop("Epsilon.mat and Xi.mat in inconsistent state")
        }

        # delta.kj
        if(!all(Delta.mat[Y.mat>0]==1)) stop("Some delta.ikj are not 1 when the subject has an observed count")
        if(!isTRUE(all.equal(Delta.sum.mat.chk, fnGetDeltaSumMat(B.mat.chk)))) stop("Delta.sum.mat in inconsistent state")
        # sigma.2.theta and sigma.2.xi
        # if(!all(sigma.2.theta.vec == sigma.2.theta)) stop("sigma.2.theta.vec inconsistent")
        # if(!all(sigma.2.xi.vec == sigma.2.xi)) stop("sigma.2.xi.vec inconsistent")

        if( (z %% SETTINGS$print.mcmc) == 0) {
            cat(format(Sys.time()), " z=", z, "\n", sep="")
        }
    })

    return(samples)
}
