// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef F_H
#define F_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::vec f(arma::vec phi,arma::vec gam, arma::mat alpha, arma::mat beta, arma::vec omega=zeros<vec>(1));

// This is the end of the header guard
#endif
