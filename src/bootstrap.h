// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef bootstrap_H
#define bootstrap_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::mat bootLoop(arma::mat X, int r, int B, double dt);
arma::mat bootLoopAggregated(arma::mat Z0, arma::mat Z1, int r, int B, double dt);
arma::mat bootstrapCpp(arma::mat X, int B, double dt);
arma::mat bootstrapCppAggregated(arma::mat Z0, arma::mat Z1, int B, double dt);

// This is the end of the header guard
#endif
