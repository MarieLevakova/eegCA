// Bootstrapping for critical values for Johansens rank test for a VECM
// Model is dX = (Pi*X[t-1]+Phi*D[t])dt + Sigma*dW[t], 3-dimensional...

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include "tools.h"
#include "johansen.h"
#include "johaTools.h"

using namespace Rcpp;
using namespace arma;
using namespace std;



// Bootstrapping loop, using the Johansen procedure for estimation of statistic-distribution (only for p=3 at this moment...)


// Bootstrapping loop for the standard setup - one trajectory
// // [[Rcpp::export]]
// arma::mat bootLoop(arma::mat X, int r, int B, double dt){
//   int N = X.n_rows-1;
//   int p = X.n_cols;
//
//   arma::mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot;
//   // Get residuals from model estimation with r=0, r=1 and r=2
//   if(r == 0){
//     Minit   = var(X, dt, false);
//     Pi      = zeros<mat>(p,p);
//   } else{
//     Minit   = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
//     alpha   = Minit.rows( 0         , r-1       );
//     beta    = Minit.rows( r         , 2*r-1     );
//     Pi      = alpha.t()*beta;
//   }
//
//   res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );
//   Psi     = Minit.rows( 2*r       , 2*r       );
//   Psi     = Psi.t();
//   S2      = Minit.rows( 2*r+1     , 2*r+p   );
//   S2      = (S2 + S2.t())/2;
//   res     = res.t();
//
//   X       = X.t();
//   Xboot   = zeros<mat>(p,N+1);
//
//   mat resBoot, Mboot;
//   mat out = zeros<mat>(B,1);
//
//   for(int b=0;b<B;b++){ // wild bootstrap
//     resBoot = (res-repmat(mean(res,1), 1, res.n_cols))%randn(p,N);
//     for(int n=1;n<N;n++){
//       Xboot.col(n) = Xboot.col(n-1) + (Pi*X.col(n-1) + Psi)*dt + chol(S2).t()*resBoot.col(n);
//     }
//
//     // Estimation where we set r = 1.
//     // Since we only extract test statistics, which are estimated in the same way for any rank,
//     // it does not really matter what rank we submit.
//     Mboot = vecm(Xboot.t(), 1, eye<mat>(p,p), eye<mat>(p,p), dt);
//
//     out(b,0) =  Mboot(2*1+p+1,r);
//   }
//   return out;
// }

// Bootstrapping loop for the aggregated setup - several distinct repeated trials
// [[Rcpp::export]]
arma::mat bootLoopAggregated(arma::mat Z0, arma::mat Z1, int r, int B, double dt,
                             int n, int n_epoch, arma::mat start, bool intercept, bool normalize){

  int N = Z0.n_rows;
  int p = Z0.n_cols;

  arma::mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot;

  // Get residuals from model estimation with r=0, r=1 and r=2
  if(r == 0){
    Minit   = varAggregated(Z0, Z1, dt, intercept, normalize);
    Pi      = zeros<mat>(p,p);
  } else{
    Minit   = vecmAggregated(Z0, Z1, r,eye<mat>(p,p), eye<mat>(p,p), dt, intercept, normalize);
    alpha   = Minit.rows( 0         , r-1       );
    beta    = Minit.rows( r         , 2*r-1     );
    Pi      = alpha.t()*beta;
  }

  res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );
  Psi     = Minit.rows( 2*r       , 2*r       );
  Psi     = Psi.t();
  S2      = Minit.rows( 2*r+1     , 2*r+p   );
  S2      = (S2 + S2.t())/2;
  res     = res.t();

  Z0       = Z0.t();
  Z1       = Z1.t();
  Xboot   = zeros<mat>(p,N+1);

  arma::mat resBoot, Mboot;

  arma::mat Z0boot = zeros<mat>(p,0);
  arma::mat Z1boot = zeros<mat>(p,0);

  arma::mat Z0boot_part   = zeros<mat>(p,n);
  arma::mat Z1boot_part   = zeros<mat>(p,n);

  arma::mat out = zeros<mat>(B,1);

  for(int b=0;b<B;b++){ // wild bootstrap
    resBoot = (res-repmat(mean(res,1), 1, res.n_cols))%randn(p,N);

    int n_iter = 0;
    for(int i=0; i<n_epoch; i++){

      // Simulate one trial
      Xboot.col(n_iter) = (start.row(i)).t(); //Insert the starting value of the trial
      for(int j=0; j<n; j++){
        n_iter = n_iter+1;
        Xboot.col(n_iter) = Xboot.col(n_iter-1) + (Pi*Xboot.col(n_iter-1) + Psi)*dt +
          chol(S2).t()*resBoot.col(n_iter-1);
        Z1boot_part.col(j) = Xboot.col(n_iter-1);
        Z0boot_part.col(j) = Xboot.col(n_iter)-Xboot.col(n_iter-1);
      }

      // Insert clones of the trial into the final dataset
      Z0boot = join_rows(Z0boot, Z0boot_part);
      Z1boot = join_rows(Z1boot, Z1boot_part);
    }

    Mboot = vecmAggregated(Z0boot.t(), Z1boot.t(), 1, eye<mat>(p,p),eye<mat>(p,p), dt, intercept, normalize); // r was originally set to 1

    out(b,0) =  Mboot(2*1+p+1,r);
  }

  return out;
}

// // Bootstrapping loop for the cloned setup - several distinct repeated trials, randomly cloned
// // [[Rcpp::export]]
// arma::mat bootLoopCloned(arma::mat Z0, arma::mat Z1, int r, int B, double dt,
//                              int n, int n_epoch, int n_aug, arma::mat start){
//
//   int N = Z0.n_rows;
//   int p = Z0.n_cols;
//
//   arma::mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot;
//
//   // Get residuals from model estimation with r=0, r=1 and r=2
//   if(r == 0){
//     Minit   = varAggregated(Z0, Z1, dt);
//     Pi      = zeros<mat>(p,p);
//     // res     = Minit.rows(p+2+1, p+2+N);
//     // Psi     = Minit.rows(0,0);
//     // Psi     = Psi.t();
//     // S2      = Minit.rows( 1     , p   );
//
//   } else{
//     Minit   = vecmAggregated(Z0, Z1, r,eye<mat>(p,p),eye<mat>(p,p), dt);
//     alpha   = Minit.rows( 0         , r-1       );
//     beta    = Minit.rows( r         , 2*r-1     );
//     Pi      = alpha.t()*beta;
//   }
//
//   res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );
//   Psi     = Minit.rows( 2*r       , 2*r       );
//   Psi     = Psi.t();
//   S2      = Minit.rows( 2*r+1     , 2*r+p   );
//   S2      = (S2 + S2.t())/2;
//   res     = res.t();
//
//   Z0       = Z0.t();
//   Z1       = Z1.t();
//   Xboot   = zeros<mat>(p,N+1);
//
//   arma::mat resBoot, Mboot, Mboot_cloned;
//
//   arma::mat Z0boot_cloned = zeros<mat>(p,0);
//   arma::mat Z1boot_cloned = zeros<mat>(p,0);
//
//   arma::mat Z0boot_part   = zeros<mat>(p,n);
//   arma::mat Z1boot_part   = zeros<mat>(p,n);
//
//   arma::mat out = zeros<mat>(B,1);
//
//   for(int b=0;b<B;b++){ // wild bootstrap
//     resBoot = (res-repmat(mean(res,1), 1, res.n_cols))%randn(p,N);
//
//     // Create a rule for cloning
//     // epoch_sample indicates which trial should be cloned and how many times
//     arma::mat epoch_sample                = randu<mat>(n_epoch*n_aug);
//     epoch_sample                          = round((n_epoch-1)*epoch_sample + 1);
//     epoch_sample.submat(0,0,n_epoch-1,0)  = linspace<vec>(1, n_epoch, n_epoch);
//     epoch_sample                          = sort(epoch_sample);
//
//     int n_iter = 0;
//     arma::mat n_clones = zeros<mat>(n_epoch,1);
//
//     for(int i=0; i<n_epoch; i++){
//       // Simulate one trial
//       Xboot.col(n_iter) = (start.row(i)).t(); //Insert the starting value of the trial
//       for(int j=0; j<n; j++){
//         n_iter = n_iter+1;
//         Xboot.col(n_iter) = Xboot.col(n_iter-1) + (Pi*Xboot.col(n_iter-1) + Psi)*dt +
//           chol(S2).t()*resBoot.col(n_iter-1);
//         Z1boot_part.col(j) = Xboot.col(n_iter-1);
//         Z0boot_part.col(j) = Xboot.col(n_iter)-Xboot.col(n_iter-1);
//       }
//
//       // Insert clones of the trial into the final dataset
//       n_clones(i,0) = accu(epoch_sample == i);
//       for(int j=0; j<n_clones(i,0); j++){
//         Z0boot_cloned = join_rows(Z0boot_cloned, Z0boot_part);
//         Z1boot_cloned = join_rows(Z1boot_cloned, Z1boot_part);
//       }
//     }
//
//     Mboot_cloned = vecmAggregated(Z0boot_cloned.t(), Z1boot_cloned.t(), 1, eye<mat>(p,p),eye<mat>(p,p), dt);
//
//     out(b,0) =  Mboot_cloned(2*1+p+1,r);
//   }
//
//   return out;
// }
//
// // Bootstrapping loop for the aggregated setup - adding noise on the top of cloned trials
// // [[Rcpp::export]]
// arma::mat bootLoopNoise(arma::mat Z0, arma::mat Z1, int r, int B, double dt,
//                              int n, int n_epoch, int n_aug,
//                              arma::mat start, double epsilon = 1){
//
//   int N = Z0.n_rows;
//   int p = Z0.n_cols;
//
//   arma::mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot;
//
//   // Get residuals from model estimation with r=0, r=1 and r=2
//   if(r == 0){
//     Minit   = varAggregated(Z0, Z1, dt);
//     Pi      = zeros<mat>(p,p);
//   } else{
//     Minit   = vecmAggregated(Z0, Z1, r,eye<mat>(p,p),eye<mat>(p,p), dt);
//     alpha   = Minit.rows( 0         , r-1       );
//     beta    = Minit.rows( r         , 2*r-1     );
//     Pi      = alpha.t()*beta;
//   }
//
//   res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );
//   Psi     = Minit.rows( 2*r       , 2*r       );
//   Psi     = Psi.t();
//   S2      = Minit.rows( 2*r+1     , 2*r+p   );
//   S2      = (S2 + S2.t())/2;
//   res     = res.t();
//
//   Z0       = Z0.t();
//   Z1       = Z1.t();
//   Xboot     = zeros<mat>(p,N+1);
//
//   arma::mat resBoot, Mboot_noise;
//
//   arma::mat Z0boot_noise = zeros<mat>(p,0);
//   arma::mat Z1boot_noise = zeros<mat>(p,0);
//
//   arma::mat Xboot_part_noise    = zeros<mat>(p,n+1);
//   arma::mat Z0boot_part         = zeros<mat>(p,n);
//   arma::mat Z1boot_part         = zeros<mat>(p,n);
//   arma::mat Z0boot_part_noise   = zeros<mat>(p,n);
//   arma::mat Z1boot_part_noise   = zeros<mat>(p,n);
//
//   arma::mat out = zeros<mat>(B,1);
//
//   for(int b=0;b<B;b++){ // wild bootstrap
//     resBoot = (res-repmat(mean(res,1), 1, res.n_cols))%randn(p,N);
//
//     int n_iter = 0;
//     for(int i=0; i<n_epoch; i++){
//
//       // Simulate one trial
//       Xboot.col(n_iter) = (start.row(i)).t(); //Insert the starting value of the trial
//       for(int j=0; j<n; j++){
//         n_iter = n_iter+1;
//         Xboot.col(n_iter) = Xboot.col(n_iter-1) + (Pi*Xboot.col(n_iter-1) + Psi)*dt +
//           chol(S2).t()*resBoot.col(n_iter-1);
//         Z1boot_part.col(j) = Xboot.col(n_iter-1);
//         Z0boot_part.col(j) = Xboot.col(n_iter)-Xboot.col(n_iter-1);
//       }
//
//       // Insert clones of the trial with noise into the final dataset
//       // (keep one copy without the added noise)
//       Z0boot_noise = join_rows(Z0boot_noise, Z0boot_part);
//       Z1boot_noise = join_rows(Z1boot_noise, Z1boot_part);
//
//       for(int j=1; j<n_aug; j++){
//         Xboot_part_noise = join_rows(Z1boot_part, Xboot.col(n_iter)) + epsilon*randn(p, n+1);
//         for(int j=0; j<n; j++){
//           Z1boot_part_noise.col(j) = Xboot_part_noise.col(j);
//           Z0boot_part_noise.col(j) = Xboot_part_noise.col(j+1) - Xboot_part_noise.col(j);
//         }
//
//         Z0boot_noise = join_rows(Z0boot_noise, Z0boot_part);
//         Z1boot_noise = join_rows(Z1boot_noise, Z1boot_part);
//       }
//     }
//
//     Mboot_noise = vecmAggregated(Z0boot_noise.t(), Z1boot_noise.t(), 1, eye<mat>(p,p),eye<mat>(p,p), dt);
//
//     out(b,0) =  Mboot_noise(2*1+p+1,r);
//   }
//
//   return out;
// }
//
// // Bootstrapping loop for the aggregated setup - randomly shifting the cloned trials
// // [[Rcpp::export]]
// arma::mat bootLoopNoise2(arma::mat Z0, arma::mat Z1, int r, int B, double dt,
//                              int n, int n_epoch, int n_aug,
//                              arma::mat start, double epsilon = 1){
//
//   int N = Z0.n_rows;
//   int p = Z0.n_cols;
//
//   arma::mat Minit, res, alpha, beta, Pi, Psi, S2, Xboot;
//
//   // Get residuals from model estimation with r=0, r=1 and r=2
//   if(r == 0){
//     Minit   = varAggregated(Z0, Z1, dt);
//     Pi      = zeros<mat>(p,p);
//
//   } else{
//     Minit   = vecmAggregated(Z0, Z1, r,eye<mat>(p,p),eye<mat>(p,p), dt);
//     alpha   = Minit.rows( 0         , r-1       );
//     beta    = Minit.rows( r         , 2*r-1     );
//     Pi      = alpha.t()*beta;
//   }
//
//   res     = Minit.rows( 2*r+p+2+1 , 2*r+p+2+N );
//   Psi     = Minit.rows( 2*r       , 2*r       );
//   Psi     = Psi.t();
//   S2      = Minit.rows( 2*r+1     , 2*r+p   );
//   S2      = (S2 + S2.t())/2;
//   res     = res.t();
//
//   Z0       = Z0.t();
//   Z1       = Z1.t();
//   Xboot   = zeros<mat>(p,N+1);
//
//   arma::mat resBoot, Mboot_noise2;
//
//   arma::mat Z0boot_noise2 = zeros<mat>(p,0);
//   arma::mat Z1boot_noise2 = zeros<mat>(p,0);
//
//   arma::mat Z0boot_part   = zeros<mat>(p,n);
//   arma::mat Z1boot_part   = zeros<mat>(p,n);
//
//   arma::mat Z0boot_part_noise   = zeros<mat>(p,n);
//   arma::mat Z1boot_part_noise   = zeros<mat>(p,n);
//
//   arma::mat out = zeros<mat>(B,1);
//
//   for(int b=0;b<B;b++){ // wild bootstrap
//     resBoot = (res-repmat(mean(res,1), 1, res.n_cols))%randn(p,N);
//
//     int n_iter = 0;
//     for(int i=0; i<n_epoch; i++){
//
//       // Simulate one trial
//       Xboot.col(n_iter) = (start.row(i)).t(); //Insert the starting value of the trial
//       for(int j=0; j<n; j++){
//         n_iter = n_iter+1;
//         Xboot.col(n_iter) = Xboot.col(n_iter-1) + (Pi*Xboot.col(n_iter-1) + Psi)*dt +
//           chol(S2).t()*resBoot.col(n_iter-1);
//         Z1boot_part.col(j) = Xboot.col(n_iter-1);
//         Z0boot_part.col(j) = Xboot.col(n_iter)-Xboot.col(n_iter-1);
//       }
//
//       // Insert clones of the trial with noise into the final dataset
//       Z0boot_noise2 = join_rows(Z0boot_noise2, Z0boot_part);
//       Z1boot_noise2 = join_rows(Z1boot_noise2, Z1boot_part);
//
//       for(int j=1; j<n_aug; j++){
//         Z0boot_noise2 = join_rows(Z0boot_noise2, Z0boot_part);
//         Z1boot_noise2 = join_rows(Z1boot_noise2, Z1boot_part + epsilon*randn()*ones<mat>(p, n));
//       }
//     }
//
//     Mboot_noise2 = vecmAggregated(Z0boot_noise2.t(), Z1boot_noise2.t(), 1, eye<mat>(p,p),eye<mat>(p,p), dt);
//
//     out(b,0) =  Mboot_noise2(2*1+p+1,r);
//
//   }
//
//   return out;
// }
//
//
// // // The function returns ...
// // [[Rcpp::export]]
// arma::mat bootstrapCpp(arma::mat X, int B, double dt){
//
//   // Include support for restricted alpha/beta matrices!
//   int N = X.n_rows-1;
//   int p = X.n_cols;
//   arma::mat est, boot, out;
//
//   // if(r>0){
//   //   est = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
//   // } else {
//   //   est = var(X, dt);
//   // }
//   if(B > 0){
//     // Run bootstrap method to determine test-statistic distributions!
//     out = zeros<mat>(B,p);
//     for(int r=0;r<p;r++){
//       boot =  bootLoop(X, r, B, dt);
//       //mat boot = zeros<mat>(B,d);
//       out.col(r) = boot; //join_cols(est,boot);
//     }
//   } else{
//     out = zeros<mat>(1,p);//est;
//   }
//   return out;
// }
//
// // [[Rcpp::export]]
// arma::mat bootstrapCppAggregated(arma::mat Z0, arma::mat Z1, int B, double dt,
//                                  int n, int n_epoch, int n_aug){
//
//   // Include support for restricted alpha/beta matrices!
//   int N = Z0.n_rows;
//   int p = Z0.n_cols;
//   arma::mat boot, out;
//
//   // if(r>0){
//   //   est = vecm(X,r,eye<mat>(p,p),eye<mat>(p,p), dt);
//   // } else {
//   //   est = var(X, dt);
//   // }
//   if(B > 0){
//     // Run bootstrap method to determine test-statistic distributions!
//     out = zeros<mat>(B,p);
//
//     for(int r=0; r<p-1; r++){
//       boot = bootLoopAggregated(Z0, Z1, r, B, dt, n, n_epoch,
//                                 zeros<mat>(n_epoch, p));
//       out.col(r) = boot; //join_cols(est,boot);
//     }
//   } else{
//     out = zeros<mat>(1,p);//est;
//   }
//   return out;
// }
//
// // [[Rcpp::export]]
//
// int myBootstrapCpp(arma::mat X, arma::vec conf_level, int B){
//   int N = X.n_rows-1;
//   int p = X.n_cols;
//
//   arma::mat jo    = johansenCpp(X, 1, eye(p,p), eye(p,p), 1);
//   arma::mat test  = jo.row(p+3);
//
//   // arma::mat boot = zeros<mat>(B,p);
//   mat upper = zeros<mat>(conf_level.n_elem,p);
//   bool insignificant = false;
//   int r = 0;
//   while(!insignificant ){
//     arma::mat boot_r = bootLoop(X, r, B, 1);
//     // boot.col(r) = boot_r;
//     upper.col(r) = quantile(boot_r, conf_level);
//     insignificant = (test(0,r)<=upper(0,r));
//     r = r+1;
//   }
//   return r-1;
// }
//
// // [[Rcpp::export]]
//
// arma::mat boot_loop(int N, int p, int r, arma::vec conf_level, int B){
//   arma::mat J = zeros<mat>(p,p);
//   for (int i=0; i<r; i++){
//     J(i, i) = -randu();
//   }
//   arma::mat P     = randn(p,p);
//   arma::mat P_inv = P.i();
//   arma::mat Pi    = P_inv * J * P;
//   arma::mat mu    = zeros<mat>(p,1);
//   arma::mat phi0  = mu;
//
//   arma::mat X = zeros<mat>(N, p);
//   X.row(0) = phi0.t();
//
//   for(int i=1;i<N;i++){
//     arma::mat dW = arma::randn(1, p);
//     X.row(i) = X.row(i-1) + X.row(i-1) * Pi.t() + mu.t() + dW;
//   }
//
//   N = X.n_rows-1;
//
//   arma::mat jo    = johansenCpp(X, 1, eye(p,p), eye(p,p), 1);
//   arma::mat test  = jo.row(p+3);
//
//   // arma::mat boot = zeros<mat>(B,p);
//   mat upper = zeros<mat>(conf_level.n_elem,p);
//   bool insignificant = false;
//   int r_iter = 0;
//   while(!insignificant ){
//     arma::mat boot_r = bootLoop(X, r_iter, B, 1);
//     // boot.col(r) = boot_r;
//     upper.col(r_iter) = quantile(boot_r, conf_level);
//     insignificant = (test(0,r_iter)<=upper(0,r_iter));
//     r_iter = r_iter+1;
//   }
//
//   arma::mat out = {{r, r_iter-1, N+1}};
//   return out;
// }
//
// // [[Rcpp::export]]
// int boot_data(arma::mat Z0, arma::mat Z1, arma::vec conf_level, int B,
//               int n, int n_epoch, arma::mat start, int r_start = 0, double dt=1){
//   int N = Z0.n_rows;
//   int p = Z0.n_cols;
//
//   arma::mat jo    = johansenCppAggregated(Z0, Z1, 1, eye(p,p), eye(p,p), dt);
//   arma::mat test  = jo.row(p+3);
//
//   // arma::mat boot = zeros<mat>(B,p);
//   mat upper = zeros<mat>(conf_level.n_elem,p);
//   mat boot = zeros<mat>(B,p);
//   bool insignificant = false;
//   int r_iter = r_start;
//   while(!insignificant & (r_iter<p)){
//     arma::mat boot_r = bootLoopAggregated(Z0, Z1, r_iter, B, dt, n, n_epoch, start);
//     boot.col(r_iter) = boot_r;
//     upper.col(r_iter) = quantile(boot_r, conf_level);
//     insignificant = (test(0,r_iter)<=upper(0,r_iter));
//     r_iter = r_iter+1;
//   }
//
//   if(!insignificant){
//     r_iter = r_iter+1;
//   }
//
//   // return boot;
//   return r_iter-1;
// }
//
// // [[Rcpp::export]]
// arma::mat boot_data_all(arma::mat Z0, arma::mat Z1, arma::vec conf_level, int B,
//                         int n, int n_epoch, arma::mat start, int r_start = 0, double dt=1){
//   int N = Z0.n_rows;
//   int p = Z0.n_cols;
//
//   arma::mat jo    = johansenCppAggregated(Z0, Z1, 1, eye(p,p), eye(p,p), dt);
//   arma::mat test  = jo.row(p+3);
//
//   cout << test << endl;
//
//   // arma::mat boot = zeros<mat>(B,p);
//   mat upper = zeros<mat>(conf_level.n_elem,p);
//   mat boot = zeros<mat>(B,p);
//   bool insignificant = false;
//   int r_iter = r_start;
//   while((r_iter<p)){
//     arma::mat boot_r = bootLoopAggregated(Z0, Z1, r_iter, B, dt, n, n_epoch, start);
//     boot.col(r_iter) = boot_r;
//     upper.col(r_iter) = quantile(boot_r, conf_level);
//     insignificant = (test(0,r_iter)<=upper(0,r_iter));
//
//     r_iter = r_iter+1;
//   }
//
//   if(!insignificant){
//     r_iter = r_iter+1;
//   }
//
//   return boot;
//   // return r_iter-1;
// }
// //
// // // [[Rcpp::export]]
// // arma::mat boot_data_quantile(arma::mat Z0, arma::mat Z1, arma::vec conf_level, int B,
// //                              int n, int n_epoch, arma::mat start, int r_start = 0, double dt){
// //   int N = Z0.n_rows;
// //   int p = Z0.n_cols;
// //
// //   arma::mat jo    = johansenCppAggregated(Z0, Z1, 1, eye(p,p), eye(p,p), dt);
// //   arma::mat test  = jo.row(p+3);
// //
// //   // arma::mat boot = zeros<mat>(B,p);
// //   mat upper = zeros<mat>(conf_level.n_elem,p);
// //   mat boot = zeros<mat>(B,p);
// //   bool insignificant = false;
// //   int r_iter = r_start;
// //   while((r_iter<p)){
// //     arma::mat boot_r = bootLoopAggregated(Z0, Z1, r_iter, B, dt, n, n_epoch, start);
// //     boot.col(r_iter) = boot_r;
// //     upper.col(r_iter) = quantile(boot_r, conf_level);
// //     insignificant = (test(0,r_iter)<=upper(0,r_iter));
// //     r_iter = r_iter+1;
// //   }
// //
// //   if(!insignificant){
// //     r_iter = r_iter+1;
// //   }
// //
// //   return upper;
// //   // return r_iter-1;
// // }
//
//
