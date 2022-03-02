// Estimates the matrix Pi with a penalized rank.
// The approach of Bunea et al. 2011, the algorithm is from Chen, Dong, Chan 2013.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>     // Cpp linear algebra library
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// More precise summing function
// [[Rcpp::depends(RcppArmadillo)]]

arma::mat penRankLoop(arma::mat Ystd, arma::mat Zstd, arma::mat Pi_ols, double lambda,
                      arma::mat U, arma::vec d, arma::mat V){

  int N = Ystd.n_rows;
  int p = Ystd.n_cols;

  vec d_restricted = d;
  for(int i=0; i<p; i++){
    if(d(i)<lambda) d_restricted(i) = 0;
  }

  mat Pi_restricted = V*diagmat(d_restricted)*diagmat(1/d)*V.t()*Pi_ols;

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

  mat Ystd = (Y - ones<mat>(N,1)*meanY); //*diagmat(1/sdY);
  mat Zstd = (Z - ones<mat>(N,1)*meanZ); //*diagmat(1/sdZ);

  // Calculate the OLS estimate and its SVD
  mat Pi_ols = Ystd.t() * Zstd * (Zstd.t()*Zstd).i();

  mat U;
  vec d;
  mat V;

  arma::svd(U, d, V, Zstd*Pi_ols.t());

  // Calculate the sequence of lambdas
  double lambda_max = max(d);
  vec lambda_seq = linspace(lambda_max, lambda_min, n_lambda);
  vec crit_value = zeros<vec>(n_lambda);

  // Choose optimal lambda
  double lambda;
  double lambda_opt;
  mat Pi_restricted;

  if(crit != 0){ //lambda not chosen by crossvalidation
    for(int i=0; i<n_lambda; i++){
      lambda = lambda_seq(i);

      Pi_restricted = penRankLoop(Ystd, Zstd, Pi_ols, lambda, U, d, V);

      cout << Pi_restricted << endl;

      int k = accu(conv_to<imat>::from(Pi_restricted!=zeros<mat>(p,p)));
      mat res = Ystd - Zstd * Pi_restricted;
      mat Omega_select = (res.t() * res)/N;
      mat Omega_inv = Omega_select.i();

      double logdet_Omega;
      double sign;

      log_det(logdet_Omega, sign, Omega_select);

      if(crit==1) { // AIC
        crit_value(i) = N*p*log(2*datum::pi) + N*logdet_Omega + 2*k + trace(res*Omega_inv*res.t());
      } else if(crit==2) { // BIC
        crit_value(i) = N*p*log(2*datum::pi) + N*logdet_Omega + k*log(N) + trace(res*Omega_inv*res.t());
      } else { // HQ
        crit_value(i) = N*p*log(2*datum::pi) + N*logdet_Omega + 2*k*log(log(N)) + trace(res*Omega_inv*res.t());
      }
    }
    lambda_opt = lambda_seq(crit_value.index_min());
  } else { // CV
    vec cv = zeros<vec>(n_lambda);

    // Run crossvalidation n_cv times
    for(int cv_run=0; cv_run<n_cv; cv_run++){
      // Divide data into 5 folds
      ivec folds = randi<ivec>(N, distr_param(1, 5));
      for(int ii=0; ii<5; ii++){
        mat Ystd_cv = Ystd.rows(find(folds!=ii));
        mat Zstd_cv = Zstd.rows(find(folds!=ii));
        mat Pi_ols_cv = Ystd_cv.t()*Zstd_cv*(Zstd_cv.t()*Zstd_cv).i();

        mat U_cv;
        vec d_cv;
        mat V_cv;
        arma::svd(U_cv, d_cv, V_cv, Zstd*Pi_ols_cv.t());

        for(int i=0; i<n_lambda; i++){
          lambda = lambda_seq(i);

          mat Pi_restricted = penRankLoop(Ystd_cv, Zstd_cv, Pi_ols, lambda, U_cv, d_cv, V_cv);
          mat res = Ystd.rows(find(folds==ii)) - Zstd.rows(find(folds==ii))*Pi_restricted;
          cv(i) = cv(i) + trace(res*res.t())/(N*p*n_cv);
        }
      }
    }

    lambda_opt = lambda_seq(cv.index_min());
    crit_value = cv;
  }

  // Fit with an optimal lambda
  Pi_restricted = penRankLoop(Ystd, Zstd, Pi_ols, lambda_opt, U, d, V);

  // Final unnormalization
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

  return mat_output;
}

