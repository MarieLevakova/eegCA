#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>

using namespace std;
using namespace arma;

// Function returns cointegration relationship vector for phi-process.
arma::vec f(arma::vec phi, arma::vec gam, arma::mat alpha, arma::mat beta, arma::vec omega=zeros<vec>(1)){
  int p = phi.size();               // Number of oscillators
  arma::vec out(p, fill::zeros);
  if(omega.size()<p) // Cast first entry of omega to p-dim vector.
    omega = omega(0)*ones<vec>(p);

  out = alpha*beta.t()*(phi-omega);       // Cointegration relation!
  out = out+gam;                          // add trend
    return out;
}
