#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' Generate random numbers from a correlation matrix
//'
//' @param X matrix
//' @param n number of random numbers to be generated per col
//'
//' @return A matrix of random numbers which replicates the correlation matrix through a singular value decomposition
//'
//' @export
//'
// [[Rcpp::export]]
arma::mat randCor(arma::mat X, int n) {

  arma::mat U;
  arma::vec S;
  arma::mat V;

  svd(U, S, V, X);

  arma::mat A = zeros(size(X));

  A.diag() = sqrt(S);

  int c(X.n_cols);

  arma::mat rnd(c, n, arma::fill::randn);

  arma::mat B(((V * A) * rnd).t());

  return B;
}
