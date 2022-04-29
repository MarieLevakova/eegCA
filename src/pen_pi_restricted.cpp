// Estimates the matrix Pi of a restricted rank r in a VECM model.
// Algorithm is from Lian, Zhao:
// Parametric and semiparametric reduced-rank regression with flexible sparsity (2015)

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>     // Cpp linear algebra library
#include <stdio.h>
#include <math.h>

using namespace Rcpp;
using namespace arma;
using namespace std;

// More precise summing function
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

double accu2(arma::mat& obj){
  long double result = 0;
  for (auto iter = obj.begin(); iter != obj.end(); ++iter)
  {
    result += *iter;
  }
  return result;
}

arma::mat penPiLoop(arma::mat Ystd, arma::mat Zstd, arma::mat Pi_init, double lambda,
                    int r, int maxiter, bool w_auto, double q, arma::mat weights){

  int N = Ystd.n_rows;
  int p = Ystd.n_cols;
  vec w = zeros<vec>(maxiter);

  mat Pi_restricted = Pi_init;

  mat PI_iter = zeros<mat>(maxiter-1, p*p);
  vec objective_iter = zeros<vec>(maxiter-1);

  for(int ii = 1; ii < maxiter; ii++){
  // 1) Soft-threshold the current iterate (penalization step)
    mat PI_free = Pi_restricted;

    for(int iii = 0; iii < p; iii++){
      for(int jjj = 0; jjj < p; jjj++){
        if (PI_free(iii,jjj) > w(ii-1)*lambda*weights(iii,jjj)) {
          PI_free(iii,jjj) = PI_free(iii,jjj) - w(ii-1)*lambda*weights(iii,jjj);
        }
        else if (PI_free(iii,jjj) < -w(ii-1)*lambda*weights(iii,jjj)) {
          PI_free(iii,jjj) = PI_free(iii,jjj) + w(ii-1)*lambda*weights(iii,jjj);
        }
        else {
          PI_free(iii,jjj) = 0;
        }
      }
    }

    // 2) Gradient step
    // Find optimal step
    if(w_auto){
      mat num1_mat = (Zstd.t()*(Ystd-Zstd*PI_free)*(Ystd-Zstd*PI_free).t()*Zstd);
      vec num1_diag = num1_mat.diag();
      long double numerator1 = accu2(num1_diag);
      mat denom_mat = (Ystd-Zstd*PI_free).t()*Zstd*Zstd.t()*Zstd*Zstd.t()*(Ystd-Zstd*PI_free);
      vec denom_diag = denom_mat.diag();
      long double denominator = accu2(denom_diag);

      mat slopesmat = (Zstd.t()*(Ystd-Zstd*PI_free));

      mat break_points_all = sort(-N*PI_free / (Zstd.t()*(Ystd-Zstd*PI_free)));
      vec break_points_vec = break_points_all.elem(find_finite(break_points_all));

      double epsilon = min(abs(diff(break_points_vec)))/2;
      vec break_points = zeros<vec>(break_points_vec.n_elem + 1);
      break_points(0) = -datum::inf;
      break_points(span(1,break_points_vec.n_elem)) = break_points_vec;
      for(int iii = 1; iii<break_points.n_elem; iii++){
        double int_up = break_points(iii);
        mat signmat = sign(((int_up-epsilon)*(Zstd.t()*(Ystd-Zstd*PI_free))/N + PI_free.t()));
        mat summat = slopesmat % signmat;
        long double numerator2 = N*accu2(summat);
        double deltat_cand = (numerator1+lambda*numerator2)/denominator;
        if(deltat_cand < int_up && deltat_cand >= break_points(iii-1)){
          w[ii] = deltat_cand;
          break;
        }
      }

      // Prevent nonsensical steps
      if(w[ii]<0){
        w[ii] = 0;
      } else if(w[ii]>1) {
        w[ii] = 1;
      }
    } else {
      w[ii] = 0.001/pow(ii,q);
    }

    mat matToProj = PI_free + w[ii]*Zstd.t()*(Ystd-Zstd*PI_free);

    // 3) Projection onto space of rank-restricted matrices

    mat svd_U;
    vec svd_d;
    mat svd_V;

    arma::svd(svd_U, svd_d, svd_V, matToProj);
    Pi_restricted = zeros<mat>(p,p);
    for(int j=0; j < r; j++){
      Pi_restricted = Pi_restricted + svd_d[j] * svd_U.col(j) * svd_V.col(j).t();
    }

    PI_iter.row(ii-1) = reshape(Pi_restricted, 1, p*p);
    objective_iter[ii-1] = accu(pow(Ystd - Zstd * Pi_restricted, 2))/(2*N) + lambda*accu(abs(Pi_restricted));

  }

  // Assemble the output
  mat out = zeros<mat>(p+maxiter, p*p+2);

  // Insert PI in the output
  out(span(0, p-1), span(0, p-1)) = Pi_restricted;

  // Insert PI_iter in the output
  out(span(p, p+maxiter-2), span(0, p*p-1)) = PI_iter;

  // Insert objective_iter in the output
  out(span(p, p+maxiter-2), p*p) = objective_iter;

  // Insert w in the output
  out(span(p, p+maxiter-1), p*p+1) = w;

  // return out;
  return out;
}

// Function to estimate Pi, lambda chosen by cross-validation method
// based on average error of 1-step forecast
// INPUT:
//   yt - multivariate time series
//   nlambda - length of lambda sequence
//   lambda_min - minimum lambda in the sequence
//   r - rank of Pi
//   maxiter - number of iterations
//   crit - numerical value representing the method to choose lambda
//   dt - timestep
//   w_auto - should weights be chosen by the optimizing procedure, or a fixed sequence
//   n_cv - number of repetitions of the crossvalidation procedure
//   weights - pxp matrix of weights for adaptive lasso
// OUTPUT:
//   mat_output - a matrix containing:
//                            Pi, mu, Omega, chosen lambda, iterations of Pi,
//                            sequence of values of the objective function, weights,
//                            full sequence of lambdas, value of the criteria to choose lambda

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]

arma::mat penPiCpp(arma::mat X, int n_lambda, double lambda_min, int r,
                   int maxiter, int crit, double dt, bool w_auto, int n_cv,
                   double q, arma::mat weights){

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

  mat Ystd = (Y - ones<mat>(N,1)*meanY)*diagmat(1/sdY);
  mat Zstd = (Z - ones<mat>(N,1)*meanZ)*diagmat(1/sdZ);

  // Calculate the sequence of lambdas
  double lambda_max = (abs(Ystd.t() * Zstd)/N).max();
  vec lambda_seq = logspace(log10(lambda_max), log10(lambda_min), n_lambda);
  vec crit_value = zeros<vec>(n_lambda);

  // Choose optimal lambda
  double lambda;
  double lambda_opt;
  mat Pi_restricted;
  mat Pi_init = zeros<mat>(p,p); //(Zstd.t()*Zstd).i()*Zstd.t()*Ystd;

  if(crit != 0){ //lambda not chosen by crossvalidation
    for(int i=0; i<n_lambda; i++){
      lambda = lambda_seq(i);

      mat pen_out = penPiLoop(Ystd, Zstd, Pi_init, lambda, r, maxiter, w_auto,
                              q, weights);
      Pi_restricted = pen_out(span(0, p-1), span(0, p-1));
      Pi_init = Pi_restricted;

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

        Pi_init = (Zstd_cv.t()*Zstd_cv).i()*Zstd_cv.t()*Ystd_cv;

        for(int i=0; i<n_lambda; i++){
          lambda = lambda_seq(i);
          mat pen_out = penPiLoop(Ystd_cv, Zstd_cv, Pi_init, lambda, r, maxiter+1,
                                  w_auto, q, weights);

          Pi_restricted = pen_out(span(0,p-1), span(0,p-1));

          Pi_init = Pi_restricted;
          mat res = Ystd.rows(find(folds==ii)) - Zstd.rows(find(folds==ii))*Pi_restricted;
          cv(i) = cv(i) + trace(res*res.t())/(N*p*n_cv);
        }
      }
    }

    lambda_opt = lambda_seq(cv.index_min());
    crit_value = cv;
  }

  // Fit with an optimal lambda
  // (Go through the whole sequence of lambdas to achieve warm starts)
  Pi_init = (Zstd.t()*Zstd).i()*Zstd.t()*Ystd;
  mat PI_iter;
  mat objective_iter;
  mat Pi_lambda = zeros<mat>(n_lambda, p*p);
  vec w = zeros<vec>(maxiter);
  // vec iter_lambda = zeros<vec>(n_lambda);
  mat Pi_final = zeros<mat>(p,p);
  for(int i=0; i<n_lambda; i++){
    lambda = lambda_seq(i);

    mat pen_out = penPiLoop(Ystd, Zstd, Pi_init, lambda, r, maxiter, w_auto, q,
                            weights);
    Pi_restricted = pen_out(span(0, p-1), span(0, p-1));
    Pi_init = Pi_restricted;
    Pi_lambda.row(n_lambda-i-1) = reshape(Pi_restricted, 1, p*p);
    // iter_lambda(n_lambda-i-1) = pen_out(0,p);
    if(lambda == lambda_opt){
      Pi_final = Pi_restricted;
      PI_iter = pen_out(span(p, p+maxiter-2), span(0, p*p-1));
      objective_iter = pen_out(span(p, p+maxiter-2), p*p);
      w = pen_out(span(p, p+maxiter-1), p*p+1);
    }
  }

  // Final unnormalization
  Pi_restricted = Pi_restricted.t();
  for(int ir=0; ir<p; ir++){
    Pi_restricted.row(ir) = Pi_restricted.row(ir)*diagmat(sdY(ir)/sdZ);
  }
  mat mu_hat = meanY.t() - Pi_restricted * meanZ.t();
  mat res = Y - ones<mat>(N, 1)*mu_hat.t() - Z*Pi_restricted.t();
  mat Omega_hat = (res.t()*res)/N;

  // Assemble the output
  mat mat_output = zeros<mat>(p+maxiter+n_lambda,p*p+2);

  mat_output(span(0, p-1), span(0, p-1)) = Pi_restricted/dt;
  mat_output(span(0, p-1), p) = mu_hat/dt;
  mat_output(span(0, p-1), span(p+1, 2*p)) = Omega_hat/dt;
  mat_output(0, 2*p+1) = lambda_opt;
  mat_output(span(p, p+maxiter-2), span(0, p*p-1)) = PI_iter;
  mat_output(span(p, p+maxiter-2), p*p) = objective_iter;
  mat_output(span(p, p+maxiter-1), p*p+1) = w;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), 0) = lambda_seq;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), 1) = crit_value;
  mat_output(span(p+maxiter, p+maxiter+n_lambda-1), span(2, p*p+1)) = Pi_lambda;

  return mat_output;
}

// // [[Rcpp::export]]
//
// arma::mat penPiCppAggregated(arma::mat Y, arma::mat Z, arma::mat Pi_init, arma::vec lambda_seq,
//                    arma::vec w, int r, float cutoff, int maxiter = 10,
//                    int crit = 0, bool lambda_default = true, bool final_zeros = true){
//
//   int n = Y.n_rows;
//   int p = Y.n_cols;
//   int cutoff_n = round(cutoff*(n+1));
//
//   // Standardize variables
//
//   mat meanY = mean(Y); // Column means of Y
//   mat meanZ = mean(Z); // Column means of Z
//
//   mat sdY = stddev(Y); // Column standard deviations of Y
//   mat sdZ = stddev(Z); // Column standard deviations of Z
//
//   mat Ystd = (Y - ones<mat>(n,1)*meanY)*diagmat(1/sdY);
//   mat Zstd = (Z - ones<mat>(n,1)*meanZ)*diagmat(1/sdZ);
//
//   mat PI_restricted = Pi_init % (sdZ.t()*(1/sdY)); // transformation of Pi into new coordinates?
//
//   // Calculate the sequence of lambdas, if it was not provided by the user
//   if(lambda_default){
//     float lambda_max = abs((Y.t()*Z)/n).max();
//     lambda_seq = linspace(0, lambda_max, 100);
//   }
//
//   int n_lambda = lambda_seq.n_elem;
//   int n_cv = Y.n_rows - cutoff_n + 2;
//
//   mat pars = zeros<mat>(p,p);
//   double lambda_opt;
//
//   if(crit==0){ // If lambda is selected by cross-validation
//
//     mat pred_sse = zeros<mat>(n_cv, n_lambda);
//     for(int i=cutoff_n-1; i<n-1; i++){
//       for(int j=0; j<n_lambda; j++){
//
//         double lambda = lambda_seq[j];
//         mat yt_train = Ystd.rows(1,i);
//         mat zt_train = Zstd.rows(1,i);
//
//         pars = penPiFixed(yt_train, zt_train, PI_restricted, lambda, w, r, maxiter,
//                           final_zeros);
//
//         pred_sse(i-cutoff_n+2,j) = sum(pow(Ystd.row(i+1)-Zstd.row(1+i)*pars.t(),2));
//       }
//     }
//     lambda_opt = lambda_seq(index_min(sum(pred_sse)));
//
//   } else { // If some of the information criteria is used
//
//     vec aic_lambdaselect = zeros<vec>(n_lambda);
//     vec bic_lambdaselect = zeros<vec>(n_lambda);
//     vec hq_lambdaselect = zeros<vec>(n_lambda);
//
//     for(int j=0; j<n_lambda; j++){
//       double lambda = lambda_seq(j);
//
//       mat pars = penPiFixed(Ystd, Zstd, PI_restricted, lambda, w, r, maxiter, final_zeros);
//
//       int k = accu(sign(abs(pars))); // count of all nonzero elements of PI
//       double loglik = 2/n * trace((Ystd - Zstd*pars.t()) * (Ystd - Zstd*pars.t()).t());
//       aic_lambdaselect(j) = loglik + 2*k;
//       bic_lambdaselect(j) = loglik + k*log(n);
//       hq_lambdaselect(j) = loglik + 2*k*log(log(n));
//     }
//
//     if(crit==1){
//       lambda_opt = lambda_seq(index_min(aic_lambdaselect));
//     } else if(crit==2){
//       lambda_opt = lambda_seq(index_min(bic_lambdaselect));
//     } else if(crit==3){
//       lambda_opt = lambda_seq(index_min(hq_lambdaselect));
//     }
//   }
//
//   // Final estimate
//   pars = penPiFixed(Ystd, Zstd, PI_restricted, lambda_opt, w, r, maxiter, final_zeros);
//
//   // Unstandardize
//   PI_restricted = PI_restricted % (sdY.t() *(1/sdZ));
//
//   mat MU_restricted = mean(Y - Z*PI_restricted.t());
//
//   mat mat_output = join_cols(PI_restricted, MU_restricted, lambda_opt*ones<mat>(1,p));
//
//   return mat_output;
//
//   // return zeros<mat>(p+2,p);
// }
