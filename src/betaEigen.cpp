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
arma::mat betaEigen(arma::mat X, int r, double dt){
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
  mat cholS;

  // Solve for the eigenvalues (and vectors)
  solveMat = inv(S11)*S10*S00_1*S01;
  // cholS = chol(B.t()*S11*B,"lower");
  // }

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  mat out = eigvec;

  return out;
}

// [[Rcpp::export]]
arma::mat betaEigenAggregated(arma::mat Z0, arma::mat Z1, int r, double dt){
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
  mat cholS;

  // Solve for the eigenvalues (and vectors)
  solveMat = inv(S11)*S10*S00_1*S01;
  // cholS = chol(B.t()*S11*B,"lower");
  // }

  cx_vec eigval_cx; // complex vector
  cx_mat eigvec_cx; // complex matrix

  eig_gen(eigval_cx, eigvec_cx, solveMat); // decomposition into eigenvalues and eigenvectors

  // C++ function returns complex vectors/matrices, so extract real parts
  vec eigval = real(eigval_cx);
  mat eigvec = real(eigvec_cx);

  // Sort by eigenvalues, descending (sort vectors first!)
  eigvec = eigvec.cols(sort_index(eigval,"descend"));
  eigval = eigval(sort_index(eigval,"descend"));

  mat out = eigvec;

  return out;
}


