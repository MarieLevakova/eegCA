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

// Auxiliary function to estimate mu with Pi and Omega fixed
arma::mat nts_muCpp(arma::mat Y, arma::mat Z, arma::mat P, arma::mat Omega){

  // Eigenvalue decomposition of Omega
  vec decomp_val;
  mat decomp_vec;
  eig_sym(decomp_val, decomp_vec, Omega);

  mat TM = decomp_vec * diagmat(1/sqrt(decomp_val)) * decomp_vec.t();
  mat TMin = decomp_vec * diagmat(sqrt(decomp_val)) * decomp_vec.t();
  mat Y_new = (Y - Z * P.t()) * TM.t(); // Creating response matrix
  mat mu = TMin * (mean(Y_new)).t();

  return mu;
}

// Auxiliary function to initialize parameters, if starting values have not been provided
arma::mat get_startingValuesCpp(arma::mat Y, arma::mat Z,
                                arma::mat alpha_init, arma::mat beta_init,
                                bool calculate_ab = true){
  int p = alpha_init.n_rows;
  int r = alpha_init.n_cols;
  int n = Y.n_rows;

  if(calculate_ab){
    mat SigmaXX = diagmat(var(Z));
    mat SigmadXdX = diagmat(var(Y));
    mat SigmaXdX = cov(Z, Y);
    mat SigmadXX = cov(Y, Z);
    mat beta_init_full = SigmaXX.i() * SigmaXdX * SigmadXdX.i() * SigmadXX;
    vec decomp_val;
    mat decomp_vec;
    eig_sym(decomp_val, decomp_vec, beta_init_full);
    beta_init = decomp_vec.cols(0,r-1);

    mat SVDinit_U;
    vec SVDinit_s;
    mat SVDinit_V;
    svd(SVDinit_U, SVDinit_s, SVDinit_V, (Z*beta_init).t() * Y);
    alpha_init = SVDinit_V * SVDinit_U.t();
  }

  mat Pi_init = alpha_init * beta_init.t();
  mat resid = Y - Z*Pi_init.t();
  mat mu_init = mean(resid);

  mat Omega_init = cov(resid);

  mat out = join_cols(alpha_init.t(), beta_init.t(), Pi_init.t(), mu_init);
  out = join_cols(out, Omega_init.t());

  return out;
}


// Penalized alpha

// [[Rcpp::export]]

arma::mat penAlphaCpp(arma::mat X, int r, arma::mat alpha_init, arma::mat beta_init,
                      int crit,
                      arma::vec rho_glasso,
                      int maxiter = 10, double conv = 0.01,
                      float cutoff = 0.8, double glmnetthresh = 1e-04,
                      bool calculate_ab = true){

  // Obtain environment containing function huge
  Rcpp::Environment huge = Environment::namespace_env("huge");
  Rcpp::Environment glmnet = Environment::namespace_env("glmnet");

  // Make functions callable from C++
  Rcpp::Function huge_r = huge["huge"];
  Rcpp::Function huge_select = huge["huge.select"];
  Rcpp::Function glmnet_r = glmnet["glmnet"];

  // Dimensions, levels and differences
  int n = X.n_rows - 1;
  int p = X.n_cols;
  int cutoff_n = round(cutoff*(n+1));

  mat Y = zeros<mat>(n,p);
  mat Z = zeros<mat>(n,p);
  for(int i=0; i<n; i++){
    Y.row(i) = X.row(i+1)-X.row(i); // Fills rows of Y with differences of yt.
    Z.row(i) = X.row(i); // Fills rows of Z with levels of yt.
  }

  // Declaration of the parameters and the output matrix
  mat MU = zeros<mat>(p,1);
  mat OMEGA = eye<mat>(p,p);
  mat out = zeros<mat>(2*r + p + 3, p);

  double lambda_opt = 0;

  if(r==0){
    int it = 0; // counter of iterations

    // Initialization of parameters
    mat mu_init = zeros<mat>(1, p);
    mat Pi_init = zeros<mat>(p, p);
    mat Omega_init = eye(p, p);

    // Obtain Mu
    mu_init = nts_muCpp(Y, Z, Pi_init, Omega_init);

    // Residuals
    mat resid = Y - ones<mat>(n,1) * mu_init;

    // Determine Omega, conditional on mu
    mat covResid = cov(resid);
    if(rho_glasso.n_elem == 1){
      Rcpp:List GLASSOfit = huge_r(Named("x", covResid), Named("lambda", rho_glasso),
                                   Named("method", "glasso"), Named("cov.output", true),
                                         Named("verbose", false));
      Rcpp::List GLASSOfitIcov = GLASSOfit["icov"];
      Omega_init = as<mat>(GLASSOfitIcov[0]);
    } else {
      Rcpp::List GLASSOfit = huge_r(Named("x", resid), Named("lambda", rho_glasso),
                                    Named("method", "glasso"), Named("cov.output", true),
                                          Named("verbose", false));
      Rcpp::List huge_BIC = huge_select(Named("est", GLASSOfit), Named("criterion", "ebic"),
                                        Named("verbose", false));
      Omega_init = as<mat>(huge_BIC["opt.icov"]);
    }

    resid = Y - ones<mat>(n,1)* mu_init.t();

    out.row(0) = mu_init.t();
    out.rows(1,p) = Omega_init.t();
    out.row(p+1) = it * ones<mat>(1,p);

  } else {

    // Initialize parameters and data matrices
    mat init = get_startingValuesCpp(Y, Z, alpha_init, beta_init, calculate_ab);

    // Initialization of the iterative procedure
    mat alpha_init = (init.rows(0,r-1)).t();
    mat beta_init = (init.rows(r,2*r-1)).t();
    mat Pi_init = (init.rows(2*r,2*r+p-1)).t();
    mat mu_init = (init.row(2*r+p)).t();
    mat Omega_init = (init.rows(2*r+p+1, 2*r+2*p)).t();

    int it = 0;
    double diff_obj = 10*conv;
    vec value_obj = ones<vec>(maxiter+1);
    mat resid = Y - ones<mat>(n,1)*(mu_init).t() - Z * beta_init * alpha_init.t();
    value_obj(0) = (1/n)*trace(resid * Omega_init * resid.t()) - log(det(Omega_init)) +
      accu(abs(alpha_init));

    // Estimate iteratively until convergence...
    while((it < maxiter) && (diff_obj > conv)){

      // Estimation of alpha

      mat Y_new = Y - ones<mat>(n,1)*mu_init.t(); // Y matrix after subtracting intercept - new response matrix
      mat Y_new_sd = stddev(Y_new);
      mat Y_new_scaled = Y_new / (ones<mat>(n,1)*Y_new_sd);
      mat Z_new = Z * beta_init; // new "design" matrix

      for(int i_r=0; i_r<p; i_r++){
        Rcpp::List determine_lambdasequence = glmnet_r(Named("y", Y_new_scaled.col(i_r)),
                                                       Named("x", Z_new),
                                                       Named("intercept", false),
                                                       Named("family", "gaussian"),
                                                       Named("thresh", 1e-4));
        mat lambda_restricted = as<vec>(determine_lambdasequence["lambda"]);
        lambda_restricted = lambda_restricted.t();
        sp_mat alpha_all = as<sp_mat>(determine_lambdasequence["beta"]);
        vec df_all = as<vec>(determine_lambdasequence["df"]);

        if(crit == 0){
          // Time series cross-validation to determine value of the tuning parameter lambda_alpha

          // Declaration of cross-validation variables

          int n_alpha = lambda_restricted.n_cols;
          mat CVscore_alpha = 100000*ones<mat>(Z_new.n_rows - cutoff_n, lambda_restricted.n_cols);


          // Loop to calculate cross-validation score
          for(int i=cutoff_n; i<Z.n_rows-1; i++){

            // Training Data
            mat Ytrain = Y_new(span(0,i-1), i_r);
            mat Ytrain_sd = stddev(Ytrain);
            mat Ytrain_scaled = Ytrain / (ones<mat>(i,1)*Ytrain_sd);
            mat Ztrain = Z_new.rows(0,i-1);

            // Test Data
            double Ytest = Y_new_scaled(i,i_r);
            mat Ztest = Z_new.row(i);

            // Estimation
            Rcpp::List ALPHA_all = glmnet_r(Named("y", Ytrain_scaled),
                                            Named("x", Ztrain),
                                            Named("lambda", lambda_restricted),
                                            Named("intercept", false),
                                            Named("family", "gaussian"),
                                            Named("thresh", 1e-4));
            sp_mat ALPHA_all_beta_sparse = as<sp_mat>(ALPHA_all["beta"]);
            mat ALPHA_all_beta = mat(ALPHA_all_beta_sparse);
            vec ALPHA_all_df = as<vec>(ALPHA_all["df"]);
            // mat ALPHA_select = ALPHA_all_beta.cols(find(ALPHA_all_df)); // ALPHA in standardized scale
            // B_BETA <- B_BETASD*Ytrain.sd

            // n_alpha = ALPHA_select.n_cols;
            // mat A_CV = 1000*ones<mat>(1, n_alpha);
            // uvec j_available = find(ALPHA_all_df);

            for(int j=0; j<n_alpha; j++){
              // A_CV(1,j) = CVscoreRegCpp(ALPHA_select.col(j), Ytest, Ztest, 1);
              // int jj = j_available(j);
              CVscore_alpha(i-cutoff_n,j) = mean(pow(Ytest-Ztest*ALPHA_all_beta.col(j), 2));
              // CVscore_alpha(i-cutoff_n,jj) = mean(pow(Ytest-Ztest*ALPHA_select.col(j), 2));
            }
            // CVscore_alpha(i-cutoff_n, find(ALPHA_all_df>0)) = A_CV;
          }

          // LogicalVector CVscore_AVAILABLE_ids(n_alpha);
          // int n_alpha_available = 0;
          // mat lambda_restricted_AVAILABLE = zeros<mat>(1, 0);
          // mat CVscore_AVAILABLE = zeros<mat>(CVscore_alpha.n_rows, 0);
          // for(int j=0; j<n_alpha; j++){
          //   if(!all(CVscore_alpha.col(j)==NA_REAL)){
          //     lambda_restricted_AVAILABLE = join_rows(lambda_restricted_AVAILABLE, lambda_restricted.col(j));
          //     CVscore_AVAILABLE = join_rows(CVscore_AVAILABLE, CVscore_alpha.col(j));
          //   }
          //   // CVscore_AVAILABLE_ids(j) = AVAILABLE_LAMBDACpp(CVscore_alpha.col(j));
          //   // if(CVscore_AVAILABLE_ids(j)){
          //   //   n_alpha_available = n_alpha_available + 1;
          //   // }
          // }

          // mat lambda_restricted_AVAILABLE = lambda_restricted(1,CVscore_AVAILABLE);

          mat meanCV  = mean(CVscore_alpha);
          // mat meanCV  = mean(CVscore_AVAILABLE);

          if(meanCV.n_cols==0){
            lambda_opt = 0;
          } else {
            lambda_opt = lambda_restricted((mean(CVscore_alpha)).index_min());
            // lambda_opt = lambda_restricted_AVAILABLE((mean(CVscore_AVAILABLE)).index_min());
          }
        }

        if(crit == 2){
          // Residual sum of squares

          mat resid = - Z_new * alpha_all;  // First step of calculating the residuals
          resid.each_col() += Y_new_scaled.col(i_r); // second step

          vec res_sq = zeros<vec>(alpha_all.n_cols);
          for(int j=0; j<alpha_all.n_cols; j++){
            res_sq(j) = accu(pow(resid.col(j), 2));
          }

          vec bic = n*log(res_sq/n) + df_all*log(n);
          int lambda_ind = (bic(find(df_all))).index_min();
          mat lambda_restricted_AVAILABLE = zeros<mat>(1, 0);

          for(int j=0; j<lambda_restricted.n_cols; j++){
            if(df_all(j)>0){
              lambda_restricted_AVAILABLE = join_rows(lambda_restricted_AVAILABLE, lambda_restricted.col(j));
            }
          }

          lambda_opt = lambda_restricted_AVAILABLE(0,lambda_ind);
        }


        Rcpp::List LASSOfinal = glmnet_r(Named("y", Y_new_scaled.col(i_r)), Named("x", Z_new),
                                         Named("intercept", false), Named("lambda", lambda_opt),
                                         Named("family", "gaussian"), Named("thresh", 1e-4));
        mat LASSOfinal_alpha = mat(as<sp_mat>(LASSOfinal["beta"])).t();
        alpha_init.row(i_r) = LASSOfinal_alpha * stddev(Y_new.col(i_r));


      }

      // Estimate beta

      // Omega.inv <- solve(OMEGA.Sparse)
      mat mat1 = Z.t() * Z;
      mat mat2 = alpha_init.t() * Omega_init * alpha_init;

      mat qr1_Q;
      mat qr1_R;
      qr(qr1_Q, qr1_R, mat1);

      mat qr2_Q;
      mat qr2_R;
      qr(qr2_Q, qr2_R, mat2);

      mat inv1 = qr1_R.i() * qr1_Q.t();
      // if (try.inverse(qr.R(qr2))){
      //           beta.tilde <- matrix(NA, nrow = dim(beta.init)[1], ncol = dim(beta.init)[2])
      // } else {
      mat inv2 = solve(qr2_R, eye<mat>(r,r)) * qr2_Q.t();
      beta_init = inv1 * Z.t() * (Y - ones<mat>(n,1)*mu_init.t()) * Omega_init * alpha_init * inv2;
      // }

      Pi_init = alpha_init * beta_init.t();

      // Estimation of mu
      mu_init = mean(Y - Z * beta_init * alpha_init.t()).t();

      // // Rcout << mu_init << endl;

      // # ers2 <- Y - matrix(rep(mu.init, N), nrow = N) - Z %*% beta.init %*% t(alpha.init)


      // Estimation of Omega
      mat r0 = Y - ones<mat>(n,1)*mean(Y); // residuals from regressing Y on the intercept
      mat r1 = Z - ones<mat>(n,1)*mean(Z); // residuals from regressing Z on the intercept
      mat ers = r0 - r1 * beta_init * alpha_init.t();
      Omega_init = ers.t() * ers /n;

      // Convergence check
      resid = Y - ones<mat>(n,1) * mu_init.t() - Z * beta_init * alpha_init.t();
      value_obj(it+1) = (1/n)*trace(resid * Omega_init * resid.t()) - log(det(Omega_init)) +
        accu(abs(alpha_init));
      diff_obj = abs(value_obj(it+1) - value_obj(it))/abs(value_obj(it));
      it = it+1;
    }

    // Assembling the output
    out.row(0) = mu_init.t();
    out.rows(1,r) = alpha_init.t();
    out.rows(r+1,2*r) = beta_init.t();
    out.rows(2*r+1,2*r+p) = Omega_init.t();
    out.row(2*r+p+1) = it * ones<mat>(1,p);
    out.row(2*r+p+2) = lambda_opt * ones<mat>(1,p);
  }
  return out;
}

// [[Rcpp::export]]

arma::mat penAlphaCppAggregated(arma::mat Y, arma::mat Z, int r,
                                arma::mat alpha_init, arma::mat beta_init,
                      int crit,
                      arma::vec rho_glasso,
                      int maxiter = 10, double conv = 0.01,
                      float cutoff = 0.8, double glmnetthresh = 1e-04,
                      bool calculate_ab = true){

  // Obtain environment containing function huge
  Rcpp::Environment huge = Environment::namespace_env("huge");
  Rcpp::Environment glmnet = Environment::namespace_env("glmnet");

  // Make functions callable from C++
  Rcpp::Function huge_r = huge["huge"];
  Rcpp::Function huge_select = huge["huge.select"];
  Rcpp::Function glmnet_r = glmnet["glmnet"];

  // Dimensions, levels and differences
  int n = Y.n_rows;
  int p = Y.n_cols;
  int cutoff_n = round(cutoff*(n+1));

  // Declaration of the parameters and the output matrix
  mat MU = zeros<mat>(p,1);
  mat OMEGA = eye<mat>(p,p);
  mat out = zeros<mat>(2*r + p + 3, p);

  double lambda_opt = 0;

  if(r==0){
    int it = 0; // counter of iterations

    // Initialization of parameters
    mat mu_init = zeros<mat>(1, p);
    mat Pi_init = zeros<mat>(p, p);
    mat Omega_init = eye(p, p);

    // Obtain Mu
    mu_init = nts_muCpp(Y, Z, Pi_init, Omega_init);

    // Residuals
    mat resid = Y - ones<mat>(n,1) * mu_init;

    // Determine Omega, conditional on mu
    mat covResid = cov(resid);
    if(rho_glasso.n_elem == 1){
      Rcpp:List GLASSOfit = huge_r(Named("x", covResid), Named("lambda", rho_glasso),
                                   Named("method", "glasso"), Named("cov.output", true),
                                         Named("verbose", false));
      Rcpp::List GLASSOfitIcov = GLASSOfit["icov"];
      Omega_init = as<mat>(GLASSOfitIcov[0]);
    } else {
      Rcpp::List GLASSOfit = huge_r(Named("x", resid), Named("lambda", rho_glasso),
                                    Named("method", "glasso"), Named("cov.output", true),
                                          Named("verbose", false));
      Rcpp::List huge_BIC = huge_select(Named("est", GLASSOfit), Named("criterion", "ebic"),
                                        Named("verbose", false));
      Omega_init = as<mat>(huge_BIC["opt.icov"]);
    }

    resid = Y - ones<mat>(n,1)* mu_init.t();

    out.row(0) = mu_init.t();
    out.rows(1,p) = Omega_init.t();
    out.row(p+1) = it * ones<mat>(1,p);

  } else {

    // Initialize parameters and data matrices
    mat init = get_startingValuesCpp(Y, Z, alpha_init, beta_init, calculate_ab);

    // Initialization of the iterative procedure
    mat alpha_init = (init.rows(0,r-1)).t();
    mat beta_init = (init.rows(r,2*r-1)).t();
    mat Pi_init = (init.rows(2*r,2*r+p-1)).t();
    mat mu_init = (init.row(2*r+p)).t();
    mat Omega_init = (init.rows(2*r+p+1, 2*r+2*p)).t();

    int it = 0;
    double diff_obj = 10*conv;
    vec value_obj = ones<vec>(maxiter+1);
    mat resid = Y - ones<mat>(n,1)*(mu_init).t() - Z * beta_init * alpha_init.t();
    value_obj(0) = (1/n)*trace(resid * Omega_init * resid.t()) - log(det(Omega_init)) +
      accu(abs(alpha_init));

    // Estimate iteratively until convergence...
    while((it < maxiter) && (diff_obj > conv)){

      // Estimation of alpha

      mat Y_new = Y - ones<mat>(n,1)*mu_init.t(); // Y matrix after subtracting intercept - new response matrix
      mat Y_new_sd = stddev(Y_new);
      mat Y_new_scaled = Y_new / (ones<mat>(n,1)*Y_new_sd);
      mat Z_new = Z * beta_init; // new "design" matrix

      for(int i_r=0; i_r<p; i_r++){
        Rcpp::List determine_lambdasequence = glmnet_r(Named("y", Y_new_scaled.col(i_r)),
                                                       Named("x", Z_new),
                                                       Named("intercept", false),
                                                       Named("family", "gaussian"),
                                                       Named("thresh", 1e-4));
        mat lambda_restricted = as<vec>(determine_lambdasequence["lambda"]);
        lambda_restricted = lambda_restricted.t();
        sp_mat alpha_all = as<sp_mat>(determine_lambdasequence["beta"]);
        vec df_all = as<vec>(determine_lambdasequence["df"]);

        if(crit == 0){
          // Time series cross-validation to determine value of the tuning parameter lambda_alpha

          // Declaration of cross-validation variables

          int n_alpha = lambda_restricted.n_cols;
          mat CVscore_alpha = 100000*ones<mat>(Z_new.n_rows - cutoff_n, lambda_restricted.n_cols);


          // Loop to calculate cross-validation score
          for(int i=cutoff_n; i<Z.n_rows-1; i++){

            // Training Data
            mat Ytrain = Y_new(span(0,i-1), i_r);
            mat Ytrain_sd = stddev(Ytrain);
            mat Ytrain_scaled = Ytrain / (ones<mat>(i,1)*Ytrain_sd);
            mat Ztrain = Z_new.rows(0,i-1);

            // Test Data
            double Ytest = Y_new_scaled(i,i_r);
            mat Ztest = Z_new.row(i);

            // Estimation
            Rcpp::List ALPHA_all = glmnet_r(Named("y", Ytrain_scaled),
                                            Named("x", Ztrain),
                                            Named("lambda", lambda_restricted),
                                            Named("intercept", false),
                                            Named("family", "gaussian"),
                                            Named("thresh", 1e-4));
            sp_mat ALPHA_all_beta_sparse = as<sp_mat>(ALPHA_all["beta"]);
            mat ALPHA_all_beta = mat(ALPHA_all_beta_sparse);
            vec ALPHA_all_df = as<vec>(ALPHA_all["df"]);
            // mat ALPHA_select = ALPHA_all_beta.cols(find(ALPHA_all_df)); // ALPHA in standardized scale
            // B_BETA <- B_BETASD*Ytrain.sd

            // n_alpha = ALPHA_select.n_cols;
            // mat A_CV = 1000*ones<mat>(1, n_alpha);
            // uvec j_available = find(ALPHA_all_df);

            for(int j=0; j<n_alpha; j++){
              // A_CV(1,j) = CVscoreRegCpp(ALPHA_select.col(j), Ytest, Ztest, 1);
              // int jj = j_available(j);
              CVscore_alpha(i-cutoff_n,j) = mean(pow(Ytest-Ztest*ALPHA_all_beta.col(j), 2));
              // CVscore_alpha(i-cutoff_n,jj) = mean(pow(Ytest-Ztest*ALPHA_select.col(j), 2));
            }
            // CVscore_alpha(i-cutoff_n, find(ALPHA_all_df>0)) = A_CV;
          }

          // LogicalVector CVscore_AVAILABLE_ids(n_alpha);
          // int n_alpha_available = 0;
          // mat lambda_restricted_AVAILABLE = zeros<mat>(1, 0);
          // mat CVscore_AVAILABLE = zeros<mat>(CVscore_alpha.n_rows, 0);
          // for(int j=0; j<n_alpha; j++){
          //   if(!all(CVscore_alpha.col(j)==NA_REAL)){
          //     lambda_restricted_AVAILABLE = join_rows(lambda_restricted_AVAILABLE, lambda_restricted.col(j));
          //     CVscore_AVAILABLE = join_rows(CVscore_AVAILABLE, CVscore_alpha.col(j));
          //   }
          //   // CVscore_AVAILABLE_ids(j) = AVAILABLE_LAMBDACpp(CVscore_alpha.col(j));
          //   // if(CVscore_AVAILABLE_ids(j)){
          //   //   n_alpha_available = n_alpha_available + 1;
          //   // }
          // }

          // mat lambda_restricted_AVAILABLE = lambda_restricted(1,CVscore_AVAILABLE);

          mat meanCV  = mean(CVscore_alpha);
          // mat meanCV  = mean(CVscore_AVAILABLE);

          if(meanCV.n_cols==0){
            lambda_opt = 0;
          } else {
            lambda_opt = lambda_restricted((mean(CVscore_alpha)).index_min());
            // lambda_opt = lambda_restricted_AVAILABLE((mean(CVscore_AVAILABLE)).index_min());
          }
        }

        if(crit == 2){
          // Residual sum of squares

          mat resid = - Z_new * alpha_all;  // First step of calculating the residuals
          resid.each_col() += Y_new_scaled.col(i_r); // second step

          vec res_sq = zeros<vec>(alpha_all.n_cols);
          for(int j=0; j<alpha_all.n_cols; j++){
            res_sq(j) = accu(pow(resid.col(j), 2));
          }

          vec bic = n*log(res_sq/n) + df_all*log(n);
          int lambda_ind = (bic(find(df_all))).index_min();
          mat lambda_restricted_AVAILABLE = zeros<mat>(1, 0);

          for(int j=0; j<lambda_restricted.n_cols; j++){
            if(df_all(j)>0){
              lambda_restricted_AVAILABLE = join_rows(lambda_restricted_AVAILABLE, lambda_restricted.col(j));
            }
          }

          lambda_opt = lambda_restricted_AVAILABLE(0,lambda_ind);
        }


        Rcpp::List LASSOfinal = glmnet_r(Named("y", Y_new_scaled.col(i_r)), Named("x", Z_new),
                                         Named("intercept", false), Named("lambda", lambda_opt),
                                         Named("family", "gaussian"), Named("thresh", 1e-4));
        mat LASSOfinal_alpha = mat(as<sp_mat>(LASSOfinal["beta"])).t();
        alpha_init.row(i_r) = LASSOfinal_alpha * stddev(Y_new.col(i_r));


      }

      // Estimate beta

      // Omega.inv <- solve(OMEGA.Sparse)
      mat mat1 = Z.t() * Z;
      mat mat2 = alpha_init.t() * Omega_init * alpha_init;

      mat qr1_Q;
      mat qr1_R;
      qr(qr1_Q, qr1_R, mat1);

      mat qr2_Q;
      mat qr2_R;
      qr(qr2_Q, qr2_R, mat2);

      mat inv1 = qr1_R.i() * qr1_Q.t();
      // if (try.inverse(qr.R(qr2))){
      //           beta.tilde <- matrix(NA, nrow = dim(beta.init)[1], ncol = dim(beta.init)[2])
      // } else {
      mat inv2 = solve(qr2_R, eye<mat>(r,r)) * qr2_Q.t();
      beta_init = inv1 * Z.t() * (Y - ones<mat>(n,1)*mu_init.t()) * Omega_init * alpha_init * inv2;
      // }

      Pi_init = alpha_init * beta_init.t();

      // Estimation of mu
      mu_init = mean(Y - Z * beta_init * alpha_init.t()).t();

      // // Rcout << mu_init << endl;

      // # ers2 <- Y - matrix(rep(mu.init, N), nrow = N) - Z %*% beta.init %*% t(alpha.init)


      // Estimation of Omega
      mat r0 = Y - ones<mat>(n,1)*mean(Y); // residuals from regressing Y on the intercept
      mat r1 = Z - ones<mat>(n,1)*mean(Z); // residuals from regressing Z on the intercept
      mat ers = r0 - r1 * beta_init * alpha_init.t();
      Omega_init = ers.t() * ers /n;

      // Convergence check
      resid = Y - ones<mat>(n,1) * mu_init.t() - Z * beta_init * alpha_init.t();
      value_obj(it+1) = (1/n)*trace(resid * Omega_init * resid.t()) - log(det(Omega_init)) +
        accu(abs(alpha_init));
      diff_obj = abs(value_obj(it+1) - value_obj(it))/abs(value_obj(it));
      it = it+1;
    }

    // Assembling the output
    out.row(0) = mu_init.t();
    out.rows(1,r) = alpha_init.t();
    out.rows(r+1,2*r) = beta_init.t();
    out.rows(2*r+1,2*r+p) = Omega_init.t();
    out.row(2*r+p+1) = it * ones<mat>(1,p);
    out.row(2*r+p+2) = lambda_opt * ones<mat>(1,p);
  }
  return out;
}
