# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

betaEigen <- function(X, r, dt) {
    .Call('_eegCA_betaEigen', PACKAGE = 'eegCA', X, r, dt)
}

betaEigenAggregated <- function(Z0, Z1, r, dt) {
    .Call('_eegCA_betaEigenAggregated', PACKAGE = 'eegCA', Z0, Z1, r, dt)
}

bootLoopAggregated <- function(Z0, Z1, r, B, dt, n, n_epoch, start, intercept, normalize) {
    .Call('_eegCA_bootLoopAggregated', PACKAGE = 'eegCA', Z0, Z1, r, B, dt, n, n_epoch, start, intercept, normalize)
}

vecm <- function(X, r, A, B, dt, normalize = TRUE) {
    .Call('_eegCA_vecm', PACKAGE = 'eegCA', X, r, A, B, dt, normalize)
}

vecmAggregated <- function(Z0, Z1, r, A, B, dt, intercept, normalize) {
    .Call('_eegCA_vecmAggregated', PACKAGE = 'eegCA', Z0, Z1, r, A, B, dt, intercept, normalize)
}

varAggregated <- function(Z0, Z1, dt, intercept, normalize) {
    .Call('_eegCA_varAggregated', PACKAGE = 'eegCA', Z0, Z1, dt, intercept, normalize)
}

johansenCpp <- function(X, r, A, B, dt, normalize) {
    .Call('_eegCA_johansenCpp', PACKAGE = 'eegCA', X, r, A, B, dt, normalize)
}

johansenCppAggregated <- function(Z0, Z1, r, A, B, dt, intercept, normalize) {
    .Call('_eegCA_johansenCppAggregated', PACKAGE = 'eegCA', Z0, Z1, r, A, B, dt, intercept, normalize)
}

vecmEigen <- function(X) {
    .Call('_eegCA_vecmEigen', PACKAGE = 'eegCA', X)
}

vecmEigenAggregated <- function(Z0, Z1) {
    .Call('_eegCA_vecmEigenAggregated', PACKAGE = 'eegCA', Z0, Z1)
}

penreg_Rcpp <- function(Y, X, lambda, beta0, control) {
    .Call('_eegCA_penreg_Rcpp', PACKAGE = 'eegCA', Y, X, lambda, beta0, control)
}

penreg_Rcpp_XY <- function(XY, XX, lambda, beta0, control) {
    .Call('_eegCA_penreg_Rcpp_XY', PACKAGE = 'eegCA', XY, XX, lambda, beta0, control)
}

oscillator <- function(N, dt, z0, alpha, beta, omega, freq, lvl, S_phi, S_gam, model) {
    .Call('_eegCA_oscillator', PACKAGE = 'eegCA', N, dt, z0, alpha, beta, omega, freq, lvl, S_phi, S_gam, model)
}

penAlphaCpp <- function(X, r, alpha_init, beta_init, crit, rho_glasso, maxiter = 10L, conv = 0.01, cutoff = 0.8, glmnetthresh = 1e-04, calculate_ab = TRUE) {
    .Call('_eegCA_penAlphaCpp', PACKAGE = 'eegCA', X, r, alpha_init, beta_init, crit, rho_glasso, maxiter, conv, cutoff, glmnetthresh, calculate_ab)
}

penAlphaCppAggregated <- function(Y, Z, r, alpha_init, beta_init, crit, rho_glasso, maxiter = 10L, conv = 0.01, cutoff = 0.8, glmnetthresh = 1e-04, calculate_ab = TRUE) {
    .Call('_eegCA_penAlphaCppAggregated', PACKAGE = 'eegCA', Y, Z, r, alpha_init, beta_init, crit, rho_glasso, maxiter, conv, cutoff, glmnetthresh, calculate_ab)
}

penAdaptNuclearCpp <- function(X, n_lambda, lambda_min, crit, dt, n_cv, w_gamma, alpha) {
    .Call('_eegCA_penAdaptNuclearCpp', PACKAGE = 'eegCA', X, n_lambda, lambda_min, crit, dt, n_cv, w_gamma, alpha)
}

penNuclearCpp <- function(X, n_lambda, lambda_min, miniter, maxiter, crit, dt, n_cv, thresh, alpha) {
    .Call('_eegCA_penNuclearCpp', PACKAGE = 'eegCA', X, n_lambda, lambda_min, miniter, maxiter, crit, dt, n_cv, thresh, alpha)
}

penRankCpp <- function(X, n_lambda, lambda_min, crit, dt, n_cv, alpha) {
    .Call('_eegCA_penRankCpp', PACKAGE = 'eegCA', X, n_lambda, lambda_min, crit, dt, n_cv, alpha)
}

accu2 <- function(obj) {
    .Call('_eegCA_accu2', PACKAGE = 'eegCA', obj)
}

penPiCpp <- function(X, n_lambda, lambda_min, r, maxiter, crit, dt, w_auto, n_cv, q, weights, lambda_max) {
    .Call('_eegCA_penPiCpp', PACKAGE = 'eegCA', X, n_lambda, lambda_min, r, maxiter, crit, dt, w_auto, n_cv, q, weights, lambda_max)
}

surr_fit_Rcpp <- function(Y, X, lambda, U0, V0, WU, WV, Xtran, control, n_cv) {
    .Call('_eegCA_surr_fit_Rcpp', PACKAGE = 'eegCA', Y, X, lambda, U0, V0, WU, WV, Xtran, control, n_cv)
}

