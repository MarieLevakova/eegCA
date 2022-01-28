// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// betaEigen
arma::mat betaEigen(arma::mat X, int r, double dt);
RcppExport SEXP _eegCA_betaEigen(SEXP XSEXP, SEXP rSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(betaEigen(X, r, dt));
    return rcpp_result_gen;
END_RCPP
}
// betaEigenAggregated
arma::mat betaEigenAggregated(arma::mat Z0, arma::mat Z1, int r, double dt);
RcppExport SEXP _eegCA_betaEigenAggregated(SEXP Z0SEXP, SEXP Z1SEXP, SEXP rSEXP, SEXP dtSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    rcpp_result_gen = Rcpp::wrap(betaEigenAggregated(Z0, Z1, r, dt));
    return rcpp_result_gen;
END_RCPP
}
// bootLoopAggregated
arma::mat bootLoopAggregated(arma::mat Z0, arma::mat Z1, int r, int B, double dt, int n, int n_epoch, arma::mat start, bool intercept, bool normalize);
RcppExport SEXP _eegCA_bootLoopAggregated(SEXP Z0SEXP, SEXP Z1SEXP, SEXP rSEXP, SEXP BSEXP, SEXP dtSEXP, SEXP nSEXP, SEXP n_epochSEXP, SEXP startSEXP, SEXP interceptSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< int >::type n_epoch(n_epochSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type start(startSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(bootLoopAggregated(Z0, Z1, r, B, dt, n, n_epoch, start, intercept, normalize));
    return rcpp_result_gen;
END_RCPP
}
// vecm
arma::mat vecm(arma::mat X, int r, arma::mat A, arma::mat B, double dt, bool normalize);
RcppExport SEXP _eegCA_vecm(SEXP XSEXP, SEXP rSEXP, SEXP ASEXP, SEXP BSEXP, SEXP dtSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(vecm(X, r, A, B, dt, normalize));
    return rcpp_result_gen;
END_RCPP
}
// vecmAggregated
arma::mat vecmAggregated(arma::mat Z0, arma::mat Z1, int r, arma::mat A, arma::mat B, double dt, bool intercept, bool normalize);
RcppExport SEXP _eegCA_vecmAggregated(SEXP Z0SEXP, SEXP Z1SEXP, SEXP rSEXP, SEXP ASEXP, SEXP BSEXP, SEXP dtSEXP, SEXP interceptSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(vecmAggregated(Z0, Z1, r, A, B, dt, intercept, normalize));
    return rcpp_result_gen;
END_RCPP
}
// varAggregated
arma::mat varAggregated(arma::mat Z0, arma::mat Z1, double dt, bool intercept, bool normalize);
RcppExport SEXP _eegCA_varAggregated(SEXP Z0SEXP, SEXP Z1SEXP, SEXP dtSEXP, SEXP interceptSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(varAggregated(Z0, Z1, dt, intercept, normalize));
    return rcpp_result_gen;
END_RCPP
}
// johansenCpp
arma::mat johansenCpp(arma::mat X, int r, arma::mat A, arma::mat B, double dt, bool normalize);
RcppExport SEXP _eegCA_johansenCpp(SEXP XSEXP, SEXP rSEXP, SEXP ASEXP, SEXP BSEXP, SEXP dtSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(johansenCpp(X, r, A, B, dt, normalize));
    return rcpp_result_gen;
END_RCPP
}
// johansenCppAggregated
arma::mat johansenCppAggregated(arma::mat Z0, arma::mat Z1, int r, arma::mat A, arma::mat B, double dt, bool intercept, bool normalize);
RcppExport SEXP _eegCA_johansenCppAggregated(SEXP Z0SEXP, SEXP Z1SEXP, SEXP rSEXP, SEXP ASEXP, SEXP BSEXP, SEXP dtSEXP, SEXP interceptSEXP, SEXP normalizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::mat >::type B(BSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< bool >::type normalize(normalizeSEXP);
    rcpp_result_gen = Rcpp::wrap(johansenCppAggregated(Z0, Z1, r, A, B, dt, intercept, normalize));
    return rcpp_result_gen;
END_RCPP
}
// vecmEigen
arma::mat vecmEigen(arma::mat X);
RcppExport SEXP _eegCA_vecmEigen(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(vecmEigen(X));
    return rcpp_result_gen;
END_RCPP
}
// vecmEigenAggregated
arma::mat vecmEigenAggregated(arma::mat Z0, arma::mat Z1);
RcppExport SEXP _eegCA_vecmEigenAggregated(SEXP Z0SEXP, SEXP Z1SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Z0(Z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z1(Z1SEXP);
    rcpp_result_gen = Rcpp::wrap(vecmEigenAggregated(Z0, Z1));
    return rcpp_result_gen;
END_RCPP
}
// oscillator
arma::mat oscillator(int N, float dt, arma::vec z0, arma::mat alpha, arma::mat beta, arma::vec omega, arma::vec freq, arma::vec lvl, arma::mat S_phi, arma::mat S_gam, const char* model);
RcppExport SEXP _eegCA_oscillator(SEXP NSEXP, SEXP dtSEXP, SEXP z0SEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP omegaSEXP, SEXP freqSEXP, SEXP lvlSEXP, SEXP S_phiSEXP, SEXP S_gamSEXP, SEXP modelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< float >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type z0(z0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type omega(omegaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type freq(freqSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lvl(lvlSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S_phi(S_phiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S_gam(S_gamSEXP);
    Rcpp::traits::input_parameter< const char* >::type model(modelSEXP);
    rcpp_result_gen = Rcpp::wrap(oscillator(N, dt, z0, alpha, beta, omega, freq, lvl, S_phi, S_gam, model));
    return rcpp_result_gen;
END_RCPP
}
// penAlphaCpp
arma::mat penAlphaCpp(arma::mat X, int r, arma::mat alpha_init, arma::mat beta_init, int crit, arma::vec rho_glasso, int maxiter, double conv, float cutoff, double glmnetthresh, bool calculate_ab);
RcppExport SEXP _eegCA_penAlphaCpp(SEXP XSEXP, SEXP rSEXP, SEXP alpha_initSEXP, SEXP beta_initSEXP, SEXP critSEXP, SEXP rho_glassoSEXP, SEXP maxiterSEXP, SEXP convSEXP, SEXP cutoffSEXP, SEXP glmnetthreshSEXP, SEXP calculate_abSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_init(alpha_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< int >::type crit(critSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rho_glasso(rho_glassoSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    Rcpp::traits::input_parameter< float >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type glmnetthresh(glmnetthreshSEXP);
    Rcpp::traits::input_parameter< bool >::type calculate_ab(calculate_abSEXP);
    rcpp_result_gen = Rcpp::wrap(penAlphaCpp(X, r, alpha_init, beta_init, crit, rho_glasso, maxiter, conv, cutoff, glmnetthresh, calculate_ab));
    return rcpp_result_gen;
END_RCPP
}
// penAlphaCppAggregated
arma::mat penAlphaCppAggregated(arma::mat Y, arma::mat Z, int r, arma::mat alpha_init, arma::mat beta_init, int crit, arma::vec rho_glasso, int maxiter, double conv, float cutoff, double glmnetthresh, bool calculate_ab);
RcppExport SEXP _eegCA_penAlphaCppAggregated(SEXP YSEXP, SEXP ZSEXP, SEXP rSEXP, SEXP alpha_initSEXP, SEXP beta_initSEXP, SEXP critSEXP, SEXP rho_glassoSEXP, SEXP maxiterSEXP, SEXP convSEXP, SEXP cutoffSEXP, SEXP glmnetthreshSEXP, SEXP calculate_abSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alpha_init(alpha_initSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type beta_init(beta_initSEXP);
    Rcpp::traits::input_parameter< int >::type crit(critSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rho_glasso(rho_glassoSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    Rcpp::traits::input_parameter< float >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< double >::type glmnetthresh(glmnetthreshSEXP);
    Rcpp::traits::input_parameter< bool >::type calculate_ab(calculate_abSEXP);
    rcpp_result_gen = Rcpp::wrap(penAlphaCppAggregated(Y, Z, r, alpha_init, beta_init, crit, rho_glasso, maxiter, conv, cutoff, glmnetthresh, calculate_ab));
    return rcpp_result_gen;
END_RCPP
}
// accu2
double accu2(arma::mat& obj);
RcppExport SEXP _eegCA_accu2(SEXP objSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type obj(objSEXP);
    rcpp_result_gen = Rcpp::wrap(accu2(obj));
    return rcpp_result_gen;
END_RCPP
}
// penPiCpp
arma::mat penPiCpp(arma::mat X, int n_lambda, double lambda_min, int r, int maxiter, int crit, double dt, bool w_auto, int n_cv);
RcppExport SEXP _eegCA_penPiCpp(SEXP XSEXP, SEXP n_lambdaSEXP, SEXP lambda_minSEXP, SEXP rSEXP, SEXP maxiterSEXP, SEXP critSEXP, SEXP dtSEXP, SEXP w_autoSEXP, SEXP n_cvSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type n_lambda(n_lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_min(lambda_minSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< int >::type crit(critSEXP);
    Rcpp::traits::input_parameter< double >::type dt(dtSEXP);
    Rcpp::traits::input_parameter< bool >::type w_auto(w_autoSEXP);
    Rcpp::traits::input_parameter< int >::type n_cv(n_cvSEXP);
    rcpp_result_gen = Rcpp::wrap(penPiCpp(X, n_lambda, lambda_min, r, maxiter, crit, dt, w_auto, n_cv));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_eegCA_betaEigen", (DL_FUNC) &_eegCA_betaEigen, 3},
    {"_eegCA_betaEigenAggregated", (DL_FUNC) &_eegCA_betaEigenAggregated, 4},
    {"_eegCA_bootLoopAggregated", (DL_FUNC) &_eegCA_bootLoopAggregated, 10},
    {"_eegCA_vecm", (DL_FUNC) &_eegCA_vecm, 6},
    {"_eegCA_vecmAggregated", (DL_FUNC) &_eegCA_vecmAggregated, 8},
    {"_eegCA_varAggregated", (DL_FUNC) &_eegCA_varAggregated, 5},
    {"_eegCA_johansenCpp", (DL_FUNC) &_eegCA_johansenCpp, 6},
    {"_eegCA_johansenCppAggregated", (DL_FUNC) &_eegCA_johansenCppAggregated, 8},
    {"_eegCA_vecmEigen", (DL_FUNC) &_eegCA_vecmEigen, 1},
    {"_eegCA_vecmEigenAggregated", (DL_FUNC) &_eegCA_vecmEigenAggregated, 2},
    {"_eegCA_oscillator", (DL_FUNC) &_eegCA_oscillator, 11},
    {"_eegCA_penAlphaCpp", (DL_FUNC) &_eegCA_penAlphaCpp, 11},
    {"_eegCA_penAlphaCppAggregated", (DL_FUNC) &_eegCA_penAlphaCppAggregated, 12},
    {"_eegCA_accu2", (DL_FUNC) &_eegCA_accu2, 1},
    {"_eegCA_penPiCpp", (DL_FUNC) &_eegCA_penPiCpp, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_eegCA(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}