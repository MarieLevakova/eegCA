// Perform Johansen procedure for alpha/beta-restricted VECM models.
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional... for now.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>     // Cpp linear algebra library
#include <stdio.h>
#include <math.h>
#include "johaTools.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::mat vecm(arma::mat X, int r, arma::mat A, arma::mat B, double dt, bool normalize=true){

  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system
  //int m = A.n_cols;     // Restrictions for alpha matrix (unused)
  int s = B.n_cols;     // Restrictions for beta matrix

  X = X.t();  // transpose   // Easier to work with column vectors

  mat Z0 = zeros<mat>(p,N);  // Matrix of zeros with p rows and N columns
  mat Z1 = zeros<mat>(p,N);
  mat Z2 = ones<mat>(1,N);   // Vector of ones.

  for(int n=0;n<N;n++){
    Z0.col(n) = X.col(n+1)-X.col(n); // Fills columns of z0 with differences of X.
    Z1.col(n) = X.col(n);            // Fills columns of z1 with levels of X.
  }


  // Find M_{ij} matrices
  mat M00,M11,M22,M01,M02,M12,M22_1;
  M00 = M(Z0,Z0);  // crossproduct divided by N, operator defined in johaTools.h
  M11 = M(Z1,Z1);
  M22 = M(Z2,Z2);
  M01 = M(Z0,Z1);
  M02 = M(Z0,Z2);
  M12 = M(Z1,Z2);

  M22_1 = inv(M22);

  // Find residuals R0,R1
  mat R0,R1;
  R0 = Z0-M02*M22_1*Z2;
  R1 = Z1-M12*M22_1*Z2;

  // Find matrices S_{ij}
  mat S00,S11,S01,S10,S00_1; // Function S does the same operation as the function M
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);

  mat solveMat;
  mat cholS;

  // For restricted A: correct residuals
  mat Ap = M_perp(A);
  mat Ab = M_bar(A);
  mat Rt0,Rt1,St00,St11,St01,St10,St00_1;

  if(Ap.size() != 1){ // A.perp is not degenerate!
    Rt0 = R0-S00*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;
    Rt1 = R1-S10*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;

    St00 = S(Rt0,Rt0);
    St11 = S(Rt1,Rt1);
    St01 = S(Rt0,Rt1);
    St10 = S(Rt1,Rt0);
    St00_1 = inv(Ab.t()*S00*Ab);

    solveMat = inv(B.t()*St11*B)*B.t()*St10*Ab*St00_1*Ab.t()*St01*B;
    cholS = chol(B.t()*St11*B,"lower");

  } else{
     // Solve for the eigenvalues (and vectors)
     solveMat = inv(B.t()*S11*B)*B.t()*S10*S00_1*S01*B;
     // cholS = chol(B.t()*S11*B,"lower");
  }

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  // Normalize eigenvectors
  for(int i=0;i<s;i++){

    if(Ap.size() != 1){ // A.perp is not degenerate!
      double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*St11*B)*eigvec.col(i)));
      eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    } else{
      double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*S11*B)*eigvec.col(i)));
      eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    }
  }

  // To use cumsum for the teststats, the eigenvalues must be sorted "in reverse", this is ascending...

  // for(int i=0;i<p;i++){
  //   if(eigval(i) > 0.9999){
  //     eigval(i) = 0.9999;
  //   }
  // }

  vec testStat = -N*cumsum(log(1-eigval(sort_index(eigval,"ascend")))); // This is the standard test of H(r) versus H(p) in an unrestricted model.

  testStat = sort(testStat,"descend");
  //uvec idx = linspace<uvec>(0, r-1, 1);

  // Normalize beta using c {p x r} = (I{r x r},0{()p-r) x r})
  mat b_hat = B*eigvec.cols(0,r-1);
  if(normalize){
    mat Ir = eye<mat>(r,r);
    mat c = join_cols(Ir,zeros<mat>(p-r,r));
    //b_hat = b_hat/b_hat(0,0);
    b_hat = b_hat*inv(c.t()*b_hat);
  }

  // Estimation of alpha and other parameters
  mat BS11B_1 = inv(b_hat.t()*S11*b_hat);

  mat a_hat = A*Ab.t()*S01*b_hat*BS11B_1;
  mat Psi_hat = (M02*M22_1-a_hat*b_hat.t()*M12*M22_1); // coefficients of lags and deterministic terms, but in fact only a constant is allowed
  mat Omega_hat = S00-S01*b_hat*BS11B_1*b_hat.t()*S10;

  // Calculate residuals
  mat Pi_hat = a_hat*b_hat.t();
  mat res = zeros<mat>(p,N);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Pi_hat*(Z1.col(n))-Psi_hat;
  }

  int outRows = a_hat.n_cols+b_hat.n_cols+1+p+2+N;
  int outCols = p;
  mat out = zeros<mat>(outRows,outCols);

  // Insert alpha estimate in output
  for(int i=0; i<a_hat.n_cols;i++){
    out.row(i) = a_hat.col(i).t()/dt;
  }

  // Insert beta estimate in output
  for(int i=0;i<b_hat.n_cols;i++){
    int j = a_hat.n_cols;
    out.row(j+i) = b_hat.col(i).t();
  }
  // Insert Psi estimate in output
  int j = a_hat.n_cols+b_hat.n_cols;
  out.row(j) = Psi_hat.t()/dt;

  // Insert Omega estimate in output
  for(int i=0;i<p;i++){
    j = a_hat.n_cols+b_hat.n_cols+1;
    out.row(j+i) = Omega_hat.col(i).t()/dt;
  }

  // Insert test statistic in output
  // Append zeros to match dimensions...
  testStat = join_cols(testStat,zeros<mat>(p-s,1));
  eigval = join_cols(eigval,zeros<mat>(p-s,1));

  j = a_hat.n_cols+b_hat.n_cols+1+p;
  out.row(j) = testStat.t();
  j = a_hat.n_cols+b_hat.n_cols+1+p+1;
  out.row(j) = eigval.t();

  // Insert residuals in output
  for(int i=0;i<N;i++){
    j = a_hat.n_cols+b_hat.n_cols+1+p+2;
    out.row(j+i) = res.col(i).t();
  }
  return out;
}

// [[Rcpp::export]]
arma::mat vecmAggregated(arma::mat Z0, arma::mat Z1, int r, arma::mat A, arma::mat B,
                         double dt, bool intercept, bool normalize){
  // Find Z0, Z1 and Z2 for input data matrix X
  int N = Z0.n_rows;   // Number of observations
  int p = Z0.n_cols;     // Dimension of the system
  //int m = A.n_cols;     // Restrictions for alpha matrix (unused)
  int s = B.n_cols;     // Restrictions for beta matrix

  Z0 = Z0.t();  // transpose   // Easier to work with column vectors
  Z1 = Z1.t();
  mat Z2 = ones<mat>(1,N);   // Vector of ones.

  mat R0,R1;
  mat M00,M11,M22,M01,M02,M12,M22_1;
  if(intercept){
    // Find matrices M_{ij}
    M00 = M(Z0,Z0);  // crossproduct divided by N, operator defined in johaTools.h
    M11 = M(Z1,Z1);
    M22 = M(Z2,Z2);
    M01 = M(Z0,Z1);
    M02 = M(Z0,Z2);
    M12 = M(Z1,Z2);

    M22_1 = inv(M22);

    // Find residuals R0,R1
    R0 = Z0-M02*M22_1*Z2;
    R1 = Z1-M12*M22_1*Z2;
  } else {
    R0 = Z0;
    R1 = Z1;
  }

  // Find matrices S_{ij}
  mat S00,S11,S01,S10,S00_1; // Function S does the same operation as the function M
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);

  mat solveMat;
  mat cholS;

  // For restricted A: correct residuals
  mat Ap = M_perp(A);
  mat Ab = M_bar(A);
  mat Rt0,Rt1,St00,St11,St01,St10,St00_1;

  if(Ap.size() != 1){ // A.perp is not degenerate!
    Rt0 = R0-S00*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;
    Rt1 = R1-S10*Ap*inv(Ap.t()*S00*Ap)*Ap.t()*R0;

    St00 = S(Rt0,Rt0);
    St11 = S(Rt1,Rt1);
    St01 = S(Rt0,Rt1);
    St10 = S(Rt1,Rt0);
    St00_1 = inv(Ab.t()*S00*Ab);

    solveMat = inv(B.t()*St11*B)*B.t()*St10*Ab*St00_1*Ab.t()*St01*B;

    cholS = chol(B.t()*St11*B,"lower");

  } else{
    // Solve for the eigenvalues (and vectors)
    solveMat = inv(B.t()*S11*B)*B.t()*S10*S00_1*S01*B;
    cholS = chol(B.t()*S11*B,"lower");
  }

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  // eig_gen(eigval_cx, eigvec_cx, S10*S00_1*S01); // decomposition into eigenvalues and eigenvectors
  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors


  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  // Normalize eigenvectors
  for(int i=0;i<s;i++){

    if(Ap.size() != 1){ // A.perp is not degenerate!
      double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*St11*B)*eigvec.col(i)));
      eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    } else{
      double nf = as_scalar(sqrt(eigvec.col(i).t()*(B.t()*S11*B)*eigvec.col(i)));
      eigvec.col(i) = eigvec.col(i)/sqrt(nf);
    }
  }


  // To use cumsum for the teststats, the eigenvalues must be sorted "in reverse", this is ascending...
  vec testStat = -N*cumsum(log(1-eigval(sort_index(eigval,"ascend")))); // This is the standard test of H(r) versus H(p) in an unrestricted model.
  vec testStat2 = -N*log(1-eigval(sort_index(eigval,"ascend"))); // This is the standard test of H(r) versus H(p) in an unrestricted model.

  testStat = sort(testStat,"descend");
  testStat2 = sort(testStat2,"descend");
  //uvec idx = linspace<uvec>(0, r-1, 1);

  // Normalize beta using c {p x r} = (I{r x r},0{()p-r) x r})
  mat b_hat = B*eigvec.cols(0,r-1);
  if(normalize){
    mat Ir = eye<mat>(r,r);
    mat c = join_cols(Ir,zeros<mat>(p-r,r));
    //b_hat = b_hat/b_hat(0,0);
    b_hat = b_hat*inv(c.t()*b_hat);
  }

  // Estimation of alpha and other parameters
  mat BS11B_1 = inv(b_hat.t()*S11*b_hat);

  mat a_hat = A*Ab.t()*S01*b_hat*BS11B_1;
  mat Psi_hat;
  if(intercept){
    Psi_hat = (M02*M22_1-a_hat*b_hat.t()*M12*M22_1); // coefficients of lags and deterministic terms, but in fact only a constant is allowed
  } else {
    Psi_hat = zeros<mat>(p,1);
  }
  // mat Omega_hat = S00-S01*b_hat*BS11B_1*b_hat.t()*S10;

  // Calculate residuals
  mat Pi_hat = a_hat*b_hat.t();
  mat res = zeros<mat>(p,N);
  mat Omega_hat = zeros<mat>(p,p);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Pi_hat*(Z1.col(n))-Psi_hat;
    Omega_hat = Omega_hat + (res.col(n)*(res.col(n)).t())/N;
  }

  int outRows = a_hat.n_cols+b_hat.n_cols+1+p+2+N+1;
  int outCols = p;
  mat out = zeros<mat>(outRows,outCols);

  // Insert alpha estimate in output
  for(int i=0; i<a_hat.n_cols;i++){
    out.row(i) = a_hat.col(i).t()/dt;
  }

  // Insert beta estimate in output
  for(int i=0;i<b_hat.n_cols;i++){
    int j = a_hat.n_cols;
    out.row(j+i) = b_hat.col(i).t();
  }
  // Insert Psi estimate in output
  int j = a_hat.n_cols+b_hat.n_cols;
  out.row(j) = Psi_hat.t()/dt;

  // Insert Omega estimate in output
  for(int i=0;i<p;i++){
    j = a_hat.n_cols+b_hat.n_cols+1;
    out.row(j+i) = Omega_hat.col(i).t()/dt;
  }

  // Insert test statistic in output
  // Append zeros to match dimensions...
  testStat = join_cols(testStat,zeros<mat>(p-s,1));
  eigval = join_cols(eigval,zeros<mat>(p-s,1));

  j = a_hat.n_cols+b_hat.n_cols+1+p;
  out.row(j) = testStat.t();
  j = a_hat.n_cols+b_hat.n_cols+1+p+1;
  out.row(j) = eigval.t();

  // Insert residuals in output
  for(int i=0;i<N;i++){
    j = a_hat.n_cols+b_hat.n_cols+1+p+2;
    out.row(j+i) = res.col(i).t();
  }

  // Insert the maximum-eigenvalue test statistic in output
  // Append zeros to match dimensions...
  testStat2 = join_cols(testStat2,zeros<mat>(p-s,1));

  j = a_hat.n_cols+b_hat.n_cols+1+p+2+N;
  out.row(j) = testStat2.t();

  return out;
}


// Estimation for a VAR model (if r=0...)
arma::mat var(arma::mat X, double dt, bool normalize){
  int N = X.n_rows-1;
  int p = X.n_cols;

  X = X.t();

  mat Z0 = zeros<mat>(p,N);
  mat Psi_hat = zeros<mat>(p,1);

  // Estimation using Least Squares
  // Formula: Psi_hat = T^-1*sum_{t=1}^T dX_t (this is 3x1 dimensional)
  for(int n=0;n<N;n++){
    Z0.col(n) = X.col(n+1)-X.col(n);
    Psi_hat += Z0.col(n);
  }
  Psi_hat = Psi_hat/N;

  // Calculate residuals and covariance estimator (see Lütkepohls book, p. 75)
  mat res = zeros<mat>(p,N);
  mat Omega = zeros<mat>(p,p);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Psi_hat;
    Omega += res.col(n)*res.col(n).t();
  }
  Omega = (Omega/N)*(N/(N-1));

  // Calculate r-hypotheses statistics
  int r_tmp = 1;
  mat tmp = vecm(X.t(), r_tmp, eye<mat>(p,p), eye<mat>(p,p), dt, normalize);
  mat test  = tmp.rows(2*r_tmp+p+1,2*r_tmp+p+1);
  mat eigs  = tmp.rows(2*r_tmp+p+2,2*r_tmp+p+2);
  mat joha = join_cols(test,eigs);

  // Model estimates
  mat est = join_cols(Psi_hat.t()/dt,Omega/dt);

  // Add test statistics for r- hypotheses
  mat est_test = join_cols(est,joha);

  // Output estimates and residuals
  mat out  = join_cols(est_test,res.t());

  // Estimation using MLE
  // Not implemented yet.

  return out;
}

// [[Rcpp::export]]
arma::mat varAggregated(arma::mat Z0, arma::mat Z1, double dt, bool intercept, bool normalize){
  int N = Z0.n_rows;
  int p = Z0.n_cols;

  Z0 = Z0.t();
  Z1 = Z1.t();

  mat Psi_hat = zeros<mat>(p,1);

  // Estimation using Least Squares
  // Formula: Psi_hat = T^-1*sum_{t=1}^T dX_t (this is 3x1 dimensional)
  if(intercept){
    for(int n=0;n<N;n++){
      Psi_hat += Z0.col(n);
    }
    Psi_hat = Psi_hat/N;
  }

  // Calculate residuals and covariance estimator (see Lütkepohls book, p. 75)
  mat res = zeros<mat>(p,N);
  mat Omega = zeros<mat>(p,p);
  for(int n=0;n<N;n++){
    res.col(n) = Z0.col(n)-Psi_hat;
    Omega += res.col(n)*res.col(n).t();
  }
  Omega = (Omega/N)*(N/(N-1));


  // Calculate r-hypotheses statistics
  int r_tmp = 1;
  mat tmp = vecmAggregated(Z0.t(),Z1.t(),r_tmp, eye<mat>(p,p), eye<mat>(p,p), dt, intercept, normalize);
  mat test  = tmp.rows(2*r_tmp+p+1,2*r_tmp+p+1);
  mat eigs  = tmp.rows(2*r_tmp+p+2,2*r_tmp+p+2);
  mat joha = join_cols(test,eigs);

  // Model estimates
  mat est = join_cols(Psi_hat.t()/dt, Omega/dt);

  // Add test statistics for r- hypotheses
  mat est_test = join_cols(est,joha);

  // Output estimates and residuals
  mat out  = join_cols(est_test,res.t());

  // Estimation using MLE
  // Not implemented yet.

  return out;
}


// Johansen estimation procedure, returns eigenvalues, eigenvectors
// [[Rcpp::export]]
arma::mat johansenCpp(arma::mat X, int r, arma::mat A, arma::mat B, double dt, bool normalize){
  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system

  mat out = zeros<mat>(N,p);
  if(r > 0){
    out = vecm(X, r, A, B, dt, normalize);
  } else{
    out = var(X, dt, normalize);
  }
  return out;
}

// [[Rcpp::export]]
arma::mat johansenCppAggregated(arma::mat Z0, arma::mat Z1, int r,
                                arma::mat A, arma::mat B, double dt,
                                bool intercept, bool normalize){
  int N = Z0.n_rows;   // Number of observations
  int p = Z0.n_cols;     // Dimension of the system

  mat out = zeros<mat>(N,p);
  if(r > 0){
    out = vecmAggregated(Z0, Z1, r, A, B, dt, intercept, normalize);
  } else {
    out = varAggregated(Z0, Z1, dt, intercept, normalize);
  }
  return out;
}

// Johansen estimation procedure, returns eigenvalues, eigenvectors
// [[Rcpp::export]]

arma::mat vecmEigen(arma::mat X){
  // Find Z0, Z1 and Z2 for input data matrix X
  int N = X.n_rows-1;   // Number of observations
  int p = X.n_cols;     // Dimension of the system

  X = X.t();  // transpose   // Easier to work with column vectors

  mat Z0 = zeros<mat>(p,N);  // Matrix of zeros with p rows and N columns
  mat Z1 = zeros<mat>(p,N);
  mat Z2 = ones<mat>(1,N);   // Vector of ones.

  for(int n=0;n<N;n++){
    Z0.col(n) = X.col(n+1)-X.col(n); // Fills columns of z0 with differences of X.
    Z1.col(n) = X.col(n);            // Fills columns of z1 with levels of X.
  }

  // Find M_{ij} matrices
  mat M00,M11,M22,M01,M02,M12,M22_1;
  M00 = M(Z0,Z0);  // crossproduct divided by N, operator defined in johaTools.h
  M11 = M(Z1,Z1);
  M22 = M(Z2,Z2);
  M01 = M(Z0,Z1);
  M02 = M(Z0,Z2);
  M12 = M(Z1,Z2);

  M22_1 = inv(M22);

  // Find residuals R0,R1
  mat R0,R1;
  R0 = Z0-M02*M22_1*Z2;
  R1 = Z1-M12*M22_1*Z2;

  // Find matrices S_{ij}
  mat S00,S11,S01,S10,S00_1; // Function S does the same operation as the function M
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);

  mat solveMat;

  // Solve for the eigenvalues (and vectors)
  solveMat = inv(S11)*S10*S00_1*S01;

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  return eigvec;
}

// [[Rcpp::export]]

arma::mat vecmEigenAggregated(arma::mat Z0, arma::mat Z1){
  // Find Z0, Z1 and Z2 for input data matrix X
  int N = Z0.n_rows;   // Number of observations
  int p = Z0.n_cols;     // Dimension of the system

  Z0 = Z0.t();  // transpose   // Easier to work with column vectors
  Z1 = Z1.t();
  mat Z2 = ones<mat>(1,N);   // Vector of ones.

  // Find M_{ij} matrices
  mat M00,M11,M22,M01,M02,M12,M22_1;
  M00 = M(Z0,Z0);  // crossproduct divided by N, operator defined in johaTools.h
  M11 = M(Z1,Z1);
  M22 = M(Z2,Z2);
  M01 = M(Z0,Z1);
  M02 = M(Z0,Z2);
  M12 = M(Z1,Z2);

  M22_1 = inv(M22);

  // Find residuals R0,R1
  mat R0,R1;
  R0 = Z0-M02*M22_1*Z2;
  R1 = Z1-M12*M22_1*Z2;

  // Find matrices S_{ij}
  mat S00,S11,S01,S10,S00_1; // Function S does the same operation as the function M
  S00 = S(R0,R0);
  S11 = S(R1,R1);
  S01 = S(R0,R1);
  S10 = S(R1,R0);
  S00_1 = inv(S00);

  mat solveMat;

  // Solve for the eigenvalues (and vectors)
  solveMat = inv(S11)*S10*S00_1*S01;

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  return eigval;
}
