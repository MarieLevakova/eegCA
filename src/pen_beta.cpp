// // Estimates the matrix Pi of a restricted rank r in a VECM model.
// // Algorithm is from Lian, Zhao:
// // Parametric and semiparametric reduced-rank regression with flexible sparsity (2015)
//
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>     // Cpp linear algebra library
// #include <stdio.h>
// #include <math.h>
//
// using namespace Rcpp;
// using namespace arma;
// using namespace std;
//
// // Auxiliary function to estimate alpha with mu, beta and Omega fixed
//
// arma::mat nts_alphaCpp(arma::mat Y, arma::mat Z, arma::mat mu, int r, arma::mat Omega, arma::mat beta){
//
//   int n = Y.n_rows;
//   arma::mat Ymatrix = Y - ones<mat>(n, 1) * mu.t();
//   arma::mat Xmatrix = Z * beta;
//
//   vec decomp_val;
//   mat decomp_vec;
//   eig_sym(decomp_val, decomp_vec, Omega);
//
//   mat TM = decomp_vec * diagmat(sqrt(decomp_val)) * decomp_vec.t();
//   mat TMin = decomp_vec * diagmat(1/sqrt(decomp_val)) * decomp_vec.t();
//
//   // Singular Value Decomposition to compute ALPHA
//   mat SingDecomp_U;
//   vec SingDecomp_s;
//   mat SingDecomp_V;
//   svd(SingDecomp_U, SingDecomp_s, SingDecomop_V, Xmatrix.t() * Ymatrix * TM);
//
//   mat alpha = TMin * SingDecomp_V * SingDecomp_U;
//
//   return alpha;
//
// }
//
// // Auxiliary function to estimate beta with mu, alpha and Omega fixed
//
// arma::mat nts_betaCpp(arma::mat Y, arma::mat Z, arma::mat mu, arma::mat Omega, arma::mat alpha,
//                       int crit,
//                       arma::rowvec,
//                       arma::vec rho_glasso, float cutoff=0.8,
//                       bool intercept = false, double glmnetthresh = 1e-04){
//
//   // Obtain environment containing function huge
//   Rcpp::Environment huge = Environment::namespace_env("huge");
//   Rcpp::Environment glmnet = Environment::namespace_env("glmnet");
//
//   // Make functions callable from C++
//   Rcpp::Function huge_r = huge["huge"];
//   Rcpp::Function huge_select = huge["huge.select"];
//   Rcpp::Function glmnet_r = glmnet["glmnet"];
//
//   int n = Y.n_rows;
//   int p = Y.n_cols;
//   int r = alpha.n_cols;
//
//   // Data matrices
//   mat Ymatrix = (Y - ones<mat>(n,1) * mu.t()) * Omega * alpha;
//   mat Xmatrix = Z;
//
//   // Store estimates of cointegrating vector
//   mat BETA_Sparse = zeros<mat>(p,r);
//
//   // Perform Lasso
//   for (int i_r=0; i_r<r; i_r++){
//     // Determine each cointegrating vector by a Lasso Regression
//
//     // Standardized Response
//     mat Ymatrixsd = Ymatrix.col(i_r)/stddev(Ymatrix.col(i_r));
//
//     // Determine lambda sequence: exclude all zero-solution
//     Rcpp::List determine_lambdasequence = glmnet_r(Named("y", Ymatrixsd),
//                                                    Named("x", Xmatrix),
//                                                    Named("standardize", false),
//                                                    Named("intercept", false),
//                                                    Named("family", "gaussian"),
//                                                    Named("thresh", glmnetthresh));
//     vec lambda_all = as<vec>(determine_lambdasequence["lambda"]);
//     vec df_all = as<vec>(determine_lambdasequence["df"]);
//     mat beta_all = as<mat>(determine_lambdasequence["beta"]);
//     // uvec df_restricted = find(df_all);
//     // vec lambda_restricted = lambda_all.elem(df_restricted);
//
//     if(crit == 2){ // BIC
//       // RSS
//       mat res = Ymatrixsd - Xmatrix * beta_all;
//       mat res_sq = sum(pow(res, 2);
//
//       vec bic = n*log(res_sq/n) + df_all*log(n);
//       int lambda_ind = bic.index_min();
//       lambda_opt = lambda_all(lambda_ind);
//     }
//
//     if(crit == 0){ // CV
//       // Time series cross-validation to determine value of the tuning parameter lambda_beta
//       int cutoff_n = round(cutoff*Ymatrix.n_rows);
//       int n_beta = lambda_all.n_elem;
//       mat CVscore_beta = 1000*ones<mat>(n-cutoff_n, lambda_all.n_cols);
//
//       for (int i=cutoff_n-1; i<n-1; i++){ // Loop to calculate cross-validation score
//         // Training Data
//         mat Ytrain = Ymatrix(span(0, i-1), i_r);
//         double Ytrain_sd = stddev(Ytrain);
//         mat Ytrain_scaled = Ytrain/Ytrain_sd;
//         mat Xtrain = Xmatrix.rows(0, i);
//
//         // Test Data
//         double Ytest = Ymatrix(i+1,i_r);
//         mat Xtest = Xmatrix.row(i+1);
//
//         // Estimation
//         Rcpp::List BETA_scaled = glmnet_r(Named("y", Ytrain_scaled),
//                                           Named("x", Xtrain),
//                                           Named("lambda", lambda_all),
//                                           Named("standardize", false),
//                                           Named("intercept", false),
//                                           Named("family", "gaussian"),
//                                           Named("thresh", glmnetthresh));
//         mat B_BETASD = as<mat>(BETA_scaled["beta"]);
//         mat B_BETA = B_BETASD*Ytrain_sd;
//
//         for(int j=0; j<n_beta; j++){
//           CVscore_beta(i-cutoff_n,j) = mean(pow(Ytest - Xtest*B_BETA.col(j), 2)/Y_sd);
//         }
//       }
//
//       int lambda_ind = mean(CVscore_beta).index_min();
//       lambda_opt = lambda_all(lambda_ind);
//     }
//     Rcpp::List LASSOfinal = glmnet_r(Named("y", Ymatrixsd),
//                                      Named("x", Xmatrix),
//                                      Named("standardize", false),
//                                      Named("intercept", false),
//                                      Named("lambda", lambda_opt),
//                                      Named("family", "gaussian"),
//                                      Named("thresh", glmnetthresh))
//     mat BETA_scaled = (as<mat>(LASSOfinal["beta"])).set_size(p, 1);
//     BETA_Sparse.col(i_r) = BETA_scaled*stddev(Ymatrix.col(i_r));
//
//   } // close loop over cointegration rank
//
//   // Determine Omega, conditional on alpha, beta and gamma
//
//   mat Resid = (Y - ones<mat>(n,1) * mu.t()) - Z * BETA_Sparse * alpha.t();
//   mat covResid = cov(Resid);
//   if(rho_glasso.n_elem == 1){
//     Rcpp::List GLASSOfit = huge_r(Named("x", covResid),
//                                   Named("lambda", rho_glasso),
//                                   Named("method", "glasso"),
//                                   Named("cov.output", true),
//                                   Named("verbose", false))
//     Rcpp::List GLASSOfitIcov = GLASSOfit["icov"];
//     OMEGA = as<mat>(GLASSOfitIcov[0]);
//   } else {
//     Rcpp::List GLASSOfit = huge_r(Named("x", Resid),
//                                   Named("lambda", rho_glasso),
//                                   Named("method", "glasso"),
//                                   Named("cov.output", true),
//                                   Named("verbose", false));
//     Rcpp::List huge_BIC = huge_select(Named("est", GLASSOfit),
//                                       Named("criterion", "ebic"),
//                                       Named("verbose", false));
//     OMEGA = as<mat>(huge_BIC["opt.icov"]);
//   }
//
//   mat out = zeros<mat>(p+r, p);
//   out.rows(0,r-1) = BETA_Sparse.t();
//   out.rows(r, p+r-1) = OMEGA;
//
//   return out;
// }
//
// // Penalized beta
//
// // [[Rcpp::export]]
//
// arma::mat penBetaCpp(arma::mat X, int r, arma::mat alpha_init, arma::mat beta_init,
//                      int crit,
//                      arma::vec rho_glasso,
//                      int maxiter = 10, double conv = 0.01,
//                      float cutoff = 0.8, double glmnetthresh = 1e-04,
//                      bool lazy_alpha = false, bool calculate_ab = true){
//
//   // Obtain environment containing function huge
//   Rcpp::Environment huge = Environment::namespace_env("huge");
//   Rcpp::Environment glmnet = Environment::namespace_env("glmnet");
//
//   // Make functions callable from C++
//   Rcpp::Function huge_r = huge["huge"];
//   Rcpp::Function huge_select = huge["huge.select"];
//   Rcpp::Function glmnet_r = glmnet["glmnet"];
//
//   // Dimensions, levels and differences
//   int n = X.n_rows - 1;
//   int p = X.n_cols;
//   int cutoff_n = round(cutoff*(n+1));
//
//   mat Y = zeros<mat>(n,p);
//   mat Z = zeros<mat>(n,p);
//   for(int i=0; i<n; i++){
//     Y.row(i) = X.row(i+1)-X.row(i); // Fills rows of Y with differences of yt.
//     Z.row(i) = X.row(i); // Fills rows of Z with levels of yt.
//   }
//
//   // Declaration of the parameters and the output matrix
//   mat MU = zeros<mat>(p,1);
//   mat OMEGA = eye<mat>(p,p);
//   mat out = zeros<mat>(2*r + p + 3, p);
//
//   double lambda_opt = 0;
//
//   if(r == 0){
//     int it = 0;
//
//     // Initialization of parameters
//     mat mu_init = zeros<mat>(1, p);
//     mat Pi_init = zeros<mat>(p, p);
//     mat Omega_init = eye(p, p);
//
//     // Obtain Mu
//     mu_init = nts_muCpp(Y, Z, Pi_init, Omega_init);
//
//     // Residuals
//     mat resid = Y - ones<mat>(n,1) * mu_init;
//
//     // Determine Omega, conditional on mu
//     mat covResid = cov(resid);
//     if(rho_glasso.n_elem == 1){
//       Rcpp:List GLASSOfit = huge_r(Named("x", covResid), Named("lambda", rho_glasso),
//                                    Named("method", "glasso"), Named("cov.output", true),
//                                          Named("verbose", false));
//       Rcpp::List GLASSOfitIcov = GLASSOfit["icov"];
//       Omega_init = as<mat>(GLASSOfitIcov[0]);
//     } else {
//       Rcpp::List GLASSOfit = huge_r(Named("x", resid), Named("lambda", rho_glasso),
//                                     Named("method", "glasso"), Named("cov.output", true),
//                                           Named("verbose", false));
//       Rcpp::List huge_BIC = huge_select(Named("est", GLASSOfit), Named("criterion", "ebic"),
//                                         Named("verbose", false));
//       Omega_init = as<mat>(huge_BIC["opt.icov"]);
//     }
//
//     resid = Y - ones<mat>(n,1)* mu_init.t();
//
//     out.row(0) = mu_init.t();
//     out.rows(1,p) = Omega_init.t();
//     out.row(p+1) = it * ones<mat>(1,p);
//
//   } else {
//
//     // Starting values
//     mat init = get_startingValuesCpp(Y, Z, alpha_init, beta_init, calculate_ab);
//
//     // Initialization of the iterative procedure
//     mat alpha_init = (init.rows(0,r-1)).t();
//     mat beta_init = (init.rows(r,2*r-1)).t();
//     mat Pi_init = (init.rows(2*r,2*r+p-1)).t();
//     mat mu_init = (init.row(2*r+p)).t();
//     mat Omega_init = (init.rows(2*r+p+1, 2*r+2*p)).t();
//
//     // Convergence parameters: initialization
//     int it = 1;
//     double diff_obj = 10*conv;
//     vec value_obj = ones<vec>(maxiter+1);
//     mat resid = Y - ones<mat>(n,1)*(mu_init).t() - Z * beta_init * alpha_init.t();
//     value_obj(0) = (1/n)*trace(resid * Omega_init * resid.t()) - log(det(Omega_init)) +
//           accu(abs(beta_init)) + accu(abs(Omega_init));
//
//     while((it < max.iter) & (diff.obj > conv)){
//       // Obtain Mu
//       mat FIT1 = nts_muCpp(Y, Z, Pi_init, Omega_init);
//
//       // Obtain Alpha
//       mat FIT2 = nts_alphaCpp(Y, Z, FIT1, r, Omega_init, beta_init);
//
//       // Obtain Beta and Omega
//       mat FIT3 = nts_betaCpp(Y, Z, FIT1, Omega_init, FIT2, crit,
//                         rho_glasso, cutoff, false, glmnetthresh);
//       mat FIT3_beta = (FIT3.rows(0, r-1)).t();
//       mat FIT3_Omega = FIT3.rows(r, p+r-1);
//
//       // Check convergence
//       mat RESID = Y - ones<mat>(n,1) * FIT1.t() - Z * FIT3_beta * FIT2.t();
//       value_obj(it) = (1/n)*trace(RESID * FIT3_Omega * RESID.t()) - log(det(FIT3_Omega)) +
//         accu(abs(FIT3_beta)) + accu(abs(FIT3_Omega));
//       diff_obj = abs(value_obj(1+it) - value_obj(it))/abs(value_obj(it));
//       beta_init = FIT3_beta;
//       alpha_init = FIT2;
//       Pi_init <- alpha_init * beta_init.t();
//       Omega_init <- FIT3_Omega;
//       it = it+1;
//     }
//
//     // Assembling the output
//     out.row(0) = mu_init.t();
//     out.rows(1,r) = alpha_init.t();
//     out.rows(r+1,2*r) = beta_init.t();
//     out.rows(2*r+1,2*r+p) = Omega_init.t();
//     out.row(2*r+p+1) = it * ones<mat>(1,p);
//     out.row(2*r+p+2) = lambda_opt * ones<mat>(1,p);
//   }
//
//   return(out)
// }
