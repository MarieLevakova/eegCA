// This is start of the header guard.  ADD_H can be any unique name.  By convention, we use the name of the header file.
#ifndef G_H
#define G_H

using namespace arma;
// This is the content of the .h file, which is where the declarations go
arma::mat g(arma::vec gam,arma::vec lvl, const char* model);

// This is the end of the header guard
#endif
