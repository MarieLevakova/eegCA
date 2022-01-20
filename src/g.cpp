#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

using namespace std;
using namespace arma;

// Function g return the deterministic part of dGamma
arma::mat g(arma::vec gam, arma::vec lvl, const char* model){
  // Reshape input phi regardless of input, so matrix multiplication makes sense
  int p = gam.size(); // Number of oscillators
  arma::vec out(p, fill::zeros);

  if(lvl.size()<p){ // Cast first entry of lvl to p-dim vector.
    lvl = lvl(0)*ones<vec>(p);
  }

  if(strcmp(model,"std")==0){
    out = zeros<vec>(p);
  } else if(strcmp(model,"win")==0){
    out = (lvl-gam)%pow(gam, 2);
  } else {
    out = zeros<vec>(p);
  }

  return out;
}
