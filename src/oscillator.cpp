/* Main function winfreeOscillator integrates a stochastically oscillating system,
 * based on linearly cointegrated phases.
 *
 * This version :  07.07.2016
 *
 * Auxillary functions used (these are found in other files):
 *  phase()       : return the phase of input of a set of (x,y) coordinates
 *  amplitude()   : return the amplitude of input of a set of (x,y) coordinates
 *  polarToXT()   : return (x,y) coordinates from polar input (gam,phi)
 *  f()           : return the deterministic trend (including the cointegration relations) for the phase process
 *  g()           : return the deterministic trend for the amplitude process
 */


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "tools.h"
#include "f.h"
#include "g.h"

using namespace std;
using namespace arma;



// [[Rcpp::export]]
arma::mat oscillator( int N,              // Number of simulations
                      float dt,           // Discretization timestep
                      arma::vec z0,       // Initial values of process
                      arma::mat alpha,    // Cointegration alpha
                      arma::mat beta,     // Cointegration beta
                      arma::vec omega,    // Cointegration omega
                      arma::vec freq,     // Intrinsic frequency
                      arma::vec lvl,      // Level for gamma process
                      arma::mat S_phi,    // Noise for phi
                      arma::mat S_gam,    // Noise for gamma
                      const char* model
                      ){

    // Idenfity the number of oscillators
        int p = z0.size()/2;

        arma::mat phi(p,N+1, fill::zeros);  // Initialize phase process
        arma::mat gam(p,N+1, fill::zeros);  // Initialize amplitude process
        arma::mat Z(2*p,N+1, fill::zeros);  // Initialize observed process

        phi.col(0)  = phase(z0);
        gam.col(0)  = amplitude(z0);
        Z.col(0)    = z0;

    // Draw noise processes for phi and gamma and scale appropriately
        arma::mat dW_phi = sqrt(dt)*randn(p,N);
        arma::mat dW_gam = sqrt(dt)*randn(p,N);

    // Main loop for iterations
        for(int n=1; n<N; n++){
          // Iterate gamma process
              gam.col(n) = gam.col(n-1)+g(gam.col(n-1), lvl, model)*dt+S_gam*dW_gam.col(n);
              if(sum(gam.col(n)) != sum(gam.col(n)) ) gam.col(n) = abs(gam.col(n));
          // Iterate phi process
              if(strcmp(model,"std")==0){
                phi.col(n) = phi.col(n-1)+f(phi.col(n-1), freq, alpha, beta, omega)*dt+S_phi*dW_phi.col(n);
              } else if (strcmp(model,"win")==0){
                phi.col(n) = phi.col(n-1)+f(phi.col(n-1), gam.col(n-1), alpha, beta, omega)*dt+S_phi*dW_phi.col(n);
              } else{
                phi.col(n) = zeros<vec>(p);
              }

          // Find observed Z-process
              Z.col(n) = polarToXY(gam.col(n), phi.col(n));
        }

    // Output is Z process + timesteps, 2*p+1 rows, N+1 columns
        vec t = regspace<vec>(0,N);
        t = dt*t;
        Z.insert_rows(0,t.t());

    return Z;
}

// THE END
