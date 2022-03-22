// Estimates the matrix Pi with a penalized rank.
// The approach of Bunea et al. 2011, the algorithm is from Chen, Dong, Chan 2013.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>     // Cpp linear algebra library
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::depends(RcppArmadillo)]]

arma::mat penRankLoop(arma::mat Ystd, arma::mat Zstd, arma::mat Pi_ols, double lambda,
                      arma::mat V){

  int N = Ystd.n_rows;
  int p = Ystd.n_cols;
  // int r = d.n_elem;

  mat W = Pi_ols*V;
  mat G = V.t();

  vec obj_k = zeros<vec>(p+1);
  obj_k(0) = pow(norm(Ystd), 2);

  for(int k=1; k<p+1; k++){
    mat Pi_k = W.cols(0,k-1) * G.rows(0,k-1);
    obj_k(k) = pow(norm(Ystd-Zstd*Pi_k), 2) + lambda*k;
  }

  int k_opt = obj_k.index_min();

  mat Pi_restricted;
  if(k_opt==0){
    Pi_restricted = zeros<mat>(p,p);
  } else {
    Pi_restricted = W.cols(0,k_opt-1) * G.rows(0,k_opt-1);
  }

  // vec d_restricted = d;
  // for(int i=0; i<r; i++){
  //   if(abs(d(i))<=lambda) d_restricted(i) = 0;
  // }
  //
  // mat Pi_restricted = Pi_ols * V * diagmat(1/d) * diagmat(d_restricted) * V.t();
  //   // V*diagmat(d_restricted)*diagmat(1/d)*V.t()*Pi_ols;

  return Pi_restricted;
}

// Function to estimate Pi, lambda chosen by cross-validation method
// based on average error of 1-step forecast
// INPUT:
//   X - multivariate time series
//   nlambda - length of lambda sequence
//   lambda_min - minimum lambda in the sequence
//   crit - numerical value representing the method to choose lambda
//   dt - timestep
//   n_cv - number of repetitions of the crossvalidation procedure
// OUTPUT:
//   mat_output - a matrix containing:
//                            Pi, mu, Omega, chosen lambda,
//                            full sequence of lambdas, value of the criteria to choose lambda

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat penRankCpp(arma::mat X, int n_lambda, double lambda_min,
                     int crit, double dt, int n_cv){

  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system
  mat Y = zeros<mat>(N,p);  // Matrix of zeros with p rows and N columns
  mat Z = zeros<mat>(N,p);
  for(int n=0;n<N;n++){
    Y.row(n) = X.row(n+1)-X.row(n); // Fills rows of Y with differences of X.
    Z.row(n) = X.row(n);            // Fills rows of Z with levels of X.
  }

  // Standardize variables
  mat meanY = mean(Y); // Column means of Y
  mat meanZ = mean(Z); // Column means of Z

  mat sdY = stddev(Y); // Column standard deviations of Y
  mat sdZ = stddev(Z); // Column standard deviations of Z

  mat Ystd = (Y - ones<mat>(N,1)*meanY);//*diagmat(1/sdY);
  mat Zstd = (Z - ones<mat>(N,1)*meanZ);//*diagmat(1/sdZ);

  // Calculate the OLS estimate and its SVD
  mat Pi_ols = pinv(Zstd.t()*Zstd) * Zstd.t() * Ystd;

  cx_vec eigval;
  cx_mat V;

  eig_gen(eigval, V, (Zstd*Pi_ols).t() * (Zstd*Pi_ols));

  mat U;
  vec d;
  mat V_svd;

  arma::svd(U, d, V_svd, Zstd*Pi_ols);
  //
  // int r = d.n_elem;

  // Calculate the sequence of lambdas
  double lambda_max = max(conv_to<vec>::from(eigval));
  vec lambda_seq = logspace(log10(lambda_max), log10(lambda_min), n_lambda);
  vec crit_value = zeros<vec>(n_lambda);
  vec aic = zeros<vec>(n_lambda);
  vec bic = zeros<vec>(n_lambda);
  vec hq = zeros<vec>(n_lambda);
  mat Pi_iter = zeros<mat>(n_lambda, p*p);

  // Choose optimal lambda
  double lambda;
  double lambda_opt;
  mat Pi_restricted;

  for(int i=0; i<n_lambda; i++){
    lambda = lambda_seq(i);

    Pi_restricted = penRankLoop(Ystd, Zstd, Pi_ols, lambda, conv_to<mat>::from(V));

    Pi_iter.row(i) = reshape(Pi_restricted, 1, p*p);

    int k = accu(conv_to<imat>::from(Pi_restricted!=zeros<mat>(p,p)));
    mat res = Ystd - Zstd * Pi_restricted;
    mat Omega_select = (res.t() * res)/N;
    mat Omega_inv = pinv(Omega_select);

    double logdet_Omega;
    double sign;

    log_det(logdet_Omega, sign, Omega_select);

    aic(i) = N*p*log(2*datum::pi) + N*logdet_Omega + 2*k + trace(res*Omega_inv*res.t());
    bic(i) = N*p*log(2*datum::pi) + N*logdet_Omega + k*log(N) + trace(res*Omega_inv*res.t());
    hq(i) = N*p*log(2*datum::pi) + N*logdet_Omega + 2*k*log(log(N)) + trace(res*Omega_inv*res.t());
  }

  if(crit==1) { // AIC
    crit_value = aic;
  } else if(crit==2) { // BIC
    crit_value = bic;
  } else { // HQ
    crit_value = hq;
  }


  if(crit==0) { // CV
    vec cv = zeros<vec>(n_lambda);

    // Run crossvalidation n_cv times
    for(int cv_run=0; cv_run<n_cv; cv_run++){
      // Divide data into 5 folds
      ivec folds = randi<ivec>(N, distr_param(1, 5));
      for(int ii=0; ii<5; ii++){
        mat Ystd_cv = Ystd.rows(find(folds!=ii));
        mat Zstd_cv = Zstd.rows(find(folds!=ii));
        mat Pi_ols_cv = pinv(Zstd_cv.t()*Zstd_cv) * Zstd_cv.t()*Ystd_cv;

        cx_vec eigval_cv;
        cx_mat V_cv;

        eig_gen(eigval_cv, V_cv, (Zstd_cv*Pi_ols_cv).t() * (Zstd_cv*Pi_ols_cv));

        // mat U_cv;
        // vec d_cv;
        // mat V_cv;
        // arma::svd(U_cv, d_cv, V_cv, Zstd*Pi_ols_cv);
        //
        // int r_cv = d_cv.n_elem;

        for(int i=0; i<n_lambda; i++){
          lambda = lambda_seq(i);

          mat Pi_restricted = penRankLoop(Ystd_cv, Zstd_cv, Pi_ols, lambda, conv_to<mat>::from(V_cv));
          mat res = Ystd.rows(find(folds==ii)) - Zstd.rows(find(folds==ii))*Pi_restricted;
          cv(i) = cv(i) + trace(res*res.t())/(N*p*n_cv);
        }
      }
    }

    crit_value = cv;
  }

  lambda_opt = lambda_seq(crit_value.index_min());

  // Fit with an optimal lambda
  Pi_restricted = penRankLoop(Ystd, Zstd, Pi_ols, lambda, conv_to<mat>::from(V));

  // Final unnormalization
  Pi_restricted = Pi_restricted.t();
  // for(int ir=0; ir<p; ir++){
  //   Pi_restricted.row(ir) = Pi_restricted.row(ir)*diagmat(sdY(ir)/sdZ);
  // }
  mat mu_hat = meanY.t() - Pi_restricted * meanZ.t();
  mat res = Y - ones<mat>(N, 1)*mu_hat.t() - Z*Pi_restricted.t();
  mat Omega_hat = (res.t()*res)/N;

  // Assemble the output
  mat mat_output = zeros<mat>(p+n_lambda,p*p+2);

  mat_output(span(0, p-1), span(0, p-1)) = Pi_restricted/dt;
  mat_output(span(0, p-1), p) = mu_hat/dt;
  mat_output(span(0, p-1), span(p+1, 2*p)) = Omega_hat/dt;
  mat_output(0, 2*p+1) = lambda_opt;
  mat_output(span(p, p+n_lambda-1), 0) = lambda_seq;
  mat_output(span(p, p+n_lambda-1), 1) = crit_value;
  mat_output(span(p, p+n_lambda-1), span(2, p*p+1)) = Pi_iter;

  return mat_output;
}

