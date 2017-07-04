#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix mat(NumericVector FG, NumericVector FL, CharacterVector RU, CharacterVector col) {

  // CharacterVector col = CharacterVector::create("GLOBAL", "ESP", "POR", "UK", "MEX", "CHI", "BRA", "ARG", "URUG", "POL", "EUR", "USA");
  int c = col.length();

  int r = RU.length();
  NumericMatrix x(r, c);

  colnames(x) = col;
  rownames(x) = RU;

  x(_,0) = FG;

  colnames(x) = col;
  rownames(x) = RU;

  for (int i = 0; i < c; i++) {
    for(int j = 0; j < r; j++) {

      x(j,i) = (RU[j] == col[i]) ? FL[j] : x(j,i);

    }
  }

  return x;

}


// [[Rcpp::export]]
List fgyfl(NumericVector x, double lim) {

  int l = x.length();
  int l2 = l / 2;
  double su = 0;

  NumericVector fg(l2);
  NumericVector fl(l2);

  for (int i = 0; i < l2; i ++) {

    fg[i] = x[i];

  }

  for (int i = l2; i < l; i++) {

    fl[i - l2] = x[i];

  }

  for (int i= 0; i < l2; i++) {

    if (pow(fg[i], 2) + pow(fl[i], 2) > lim) {

      su = sqrt(pow(fg[i], 2) + pow(fl[i], 2));

      fg[i] = sqrt(lim) * fg[i] / su;
      fl[i] = sqrt(lim) * fl[i] / su;

    }
  }

  return List::create(
    _["FG"] = fg,
    _["FL"] = fl
  );

}

