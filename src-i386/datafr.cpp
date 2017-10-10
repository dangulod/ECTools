#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix mat(NumericVector FG, NumericMatrix FL, DataFrame RU, CharacterVector col) {

  // number of columns of the matrix x
  int c = col.length();

  // number of rows of the matrix x
  int r = RU.nrows();
  int cu = RU.ncol();

  // create the output matrix
  NumericMatrix x(r, c);

  // set the col & row names
  colnames(x) = col;
  rownames(x) = (CharacterVector)RU[0];

  // set the Global Factor
  x(_, 0) = FG;

  for (int k = 1; k < cu; k++) {        // Loop through number of columns of RU

    CharacterVector ru = RU[k];

    for (int i = 0; i < r; i++) {       // Loop through number of columns of  RU
      for (int j = 1; j < c; j++) {     // Loop through number of elements of col
        if (ru[i] == col[j]) x(i,j) =  FL(i, k - 1);
      }
    }
  }

  return x;

}



// [[Rcpp::export]]
NumericMatrix fgyfl(NumericVector x, double lim, CharacterVector n) {

  // length of x vector
  int l = x.length();

  // number of factors
  int ln = n.length() + 1;

  // number of elements of every factor
  int l2 = l / ln;

  CharacterVector names(ln);
  names[0] = "FG";

  for (int i = 1; i < ln; i++) {
    names[i] = n[i -1];
  }

  //
  NumericMatrix fgyfl(l2, ln);

  colnames(fgyfl) = names;

  int r = 0;
  int c = 0;

  if (l2 * ln != l)
    throw std::range_error("'x' does not fit on 'n' length");

  for (int i = 0; i < l; i++) {
    c = i % l;
    r = i - (c * l2);
    fgyfl(r,c) = x[i];
  }

  double su = 0;

  for (int i = 0; i < l2; i++) {            // loop through the matrix rows
    su = 0;
    for (int j = 0; j < ln; j++) {          // loop through the matrix columns
      su += pow(fgyfl(i,j), 2);
    }

    if (su > lim) {
      for (int j = 0; j < ln; j++) {
        fgyfl(i,j) = sqrt(lim) * fgyfl(i,j) / sqrt(su);
      }
    }
  }

  return fgyfl;

}
