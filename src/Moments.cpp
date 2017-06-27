# include <Rcpp.h>
# include <iostream>

//' The first raw moment
//'
//' @param loss vector with input
//' @param scale scale parameter
//'
//' @return The first raw moment is the mean
//'
//' @references https://en.wikipedia.org/wiki/Moment_(mathematics)
//'
//' @export
//'
// [[Rcpp::export]]
double k1(std::vector<double> x) {

  double s = 0;
  int l = x.size();
  for (int i = 0; i < l; i++) {

    s += x[i] / l;

  };

  return s;
}

//' The second central moment
//'
//' @param loss vector with input
//' @param scale scale parameter
//'
//' @return The second central moment is the variance. Its positive square root is the standard deviation σ.
//'
//' @references https://en.wikipedia.org/wiki/Moment_(mathematics)
//'
//' @export
//'
// [[Rcpp::export]]
double k2(std::vector<double> x) {

  double s = 0;
  int l = x.size();
  double mean = k1(x);

  for (int i = 0; i < l; i++) {

    s += std::pow((x[i] - mean), 2) / l;

  };

  return s;
}

//' The third central moment
//'
//' @param loss vector with input
//' @param scale scale parameter
//'
//' @return The third central moment is the measure of the lopsidedness of the distribution; any symmetric distribution will have a third
//' central moment, if defined, of zero. The normalised third central moment is called the skewness, often Gamma. A distribution that is
//' skewed to the left (the tail of the distribution is longer on the left) will have a negative skewness. A distribution that is skewed
//' to the right (the tail of the distribution is longer on the right), will have a positive skewness.
//' For distributions that are not too different from the normal distribution, the median will be somewhere near MU − Gamma * Sigma / 6;
//' the mode about MU − Gamma * Sigma / 2.
//'
//' @references https://en.wikipedia.org/wiki/Moment_(mathematics)
//'
//' @export
//'
// [[Rcpp::export]]
double k3(std::vector<double> x) {

  double s = 0;
  int l = x.size();
  double media = k1(x);

  for (int i = 0; i < l; i++) {
    s += std::pow((x[i] - media), 3);

  };

  s = s / l;

  return s;
}

//' The fourth central moment
//'
//' @param loss vector with input
//' @param scale scale parameter
//'
//' @return The fourth central moment is a measure of the heaviness of the tail of the distribution, compared to the normal distribution
//' of the same variance. Since it is the expectation of a fourth power, the fourth central moment, where defined, is always positive;
//' and except for a point distribution, it is always strictly positive. The fourth central moment of a normal distribution is 3 * Sigma 4.
//'
//' The kurtosis K is defined to be the normalised fourth central moment minus 3 (Equivalently, as in the next section, it is the fourth
//' cumulant divided by the square of the variance). Some authorities do not subtract three, but it is usually more convenient to have the
//' normal distribution at the origin of coordinates.[4][5] If a distribution has heavy tails, the kurtosis will be high (sometimes called
//' leptokurtic); conversely, light-tailed distributions (for example, bounded distributions such as the uniform) have low kurtosis (sometimes
//' called platykurtic).
//'
//' The kurtosis can be positive without limit, but K must be greater than or equal to Gamma * 2 − 2; equality only holds for binary distributions.
//' For unbounded skew distributions not too far from normal, K tends to be somewhere in the area of Gamma * 2 and 2 * Gamma * 2.
//'
//' The inequality can be proven by considering
//'
//' E[(T^2-aT-1)^2
//'
//' where T = (X − μ)/σ. This is the expectation of a square, so it is non-negative for all a; however it is also a quadratic polynomial in a.
//' Its discriminant must be non-positive, which gives the required relationship.
//'
//' @references https://en.wikipedia.org/wiki/Moment_(mathematics)
//'
//' @export
//'
// [[Rcpp::export]]
double k4(std::vector<double> x) {

  double s = 0;
  int l = x.size();
  double media = k1(x);
  double var = k2(x);

  for (int i = 0; i < l; i++) {
    s += std::pow((x[i] - media), 4) / l;

  };

  s = s - 3 * pow(var, 2);

  return s;
}

//' The fifth central moment
//'
//' @param loss vector with input
//' @param scale scale parameter
//'
//' @return High-order moments are moments beyond 4th-order moments. As with variance, skewness, and kurtosis, these are higher-order statistics,
//' involving non-linear combinations of the data, and can be used for description or estimation of further shape parameters. The higher the moment,
//' the harder it is to estimate, in the sense that larger samples are required in order to obtain estimates of similar quality. This is due to the
//' excess degrees of freedom consumed by the higher orders. Further, they can be subtle to interpret, often being most easily understood in terms of
//' lower order moments – compare the higher derivatives of jerk and jounce in physics. For example, just as the 4th-order moment (kurtosis) can be
//' interpreted as "relative importance of tails versus shoulders in causing dispersion" (for a given dispersion, high kurtosis corresponds to heavy
//' tails, while low kurtosis corresponds to broad shoulders), the 5th-order moment can be interpreted as measuring "relative importance of tails
//' versus center (mode, shoulders) in causing skew" (for a given skew, high 5th moment corresponds to heavy tail and little movement of mode, while
//' low 5th moment corresponds to more change in shoulders).
//'
//' @references https://en.wikipedia.org/wiki/Moment_(mathematics)
//'
//' @export
//'
// [[Rcpp::export]]
double k5(std::vector<double> x) {

  double s = 0;
  int l = x.size();
  double media = k1(x);
  double var = k2(x);
  double sim = k3(x);

  for (int i = 0; i < l; i++) {
    s += std::pow((x[i] - media), 5) / l;

  };

  s = s - 10 * var * sim;

  return s;
}

//' The sixth central moment
//'
//' @param loss vector with input
//' @param scale scale parameter
//'
//' @return High-order moments are moments beyond 4th-order moments. As with variance, skewness, and kurtosis, these are higher-order statistics,
//' involving non-linear combinations of the data, and can be used for description or estimation of further shape parameters. The higher the moment,
//' the harder it is to estimate, in the sense that larger samples are required in order to obtain estimates of similar quality. This is due to the
//' excess degrees of freedom consumed by the higher orders. Further, they can be subtle to interpret, often being most easily understood in terms
//' of lower order moments – compare the higher derivatives of jerk and jounce in physics. For example, just as the 4th-order moment (kurtosis) can
//' be interpreted as "relative importance of tails versus shoulders in causing dispersion" (for a given dispersion, high kurtosis corresponds to
//' heavy tails, while low kurtosis corresponds to broad shoulders), the 5th-order moment can be interpreted as measuring "relative importance of
//' tails versus center (mode, shoulders) in causing skew" (for a given skew, high 5th moment corresponds to heavy tail and little movement of mode,
//' while low 5th moment corresponds to more change in shoulders).
//'
//' @references https://en.wikipedia.org/wiki/Moment_(mathematics)
//'
//' @export
//'
// [[Rcpp::export]]
double k6(std::vector<double> x) {

  double s = 0;
  int l = x.size();
  double media = k1(x);
  double var = k2(x);
  double sim = k3(x);
  double m4 = k4(x) + 3 * pow(var, 2);

  for (int i = 0; i < l; i++) {
    s += std::pow((x[i] - media), 6) / l;

  };

  s = s - 15 * var * m4 - 10 * pow(sim, 2) + 30 * pow(var, 3);

  return s;
}

//' @export
// [[Rcpp::export]]
std::vector<double> st_cumulants(double location, double escala, double shape, double df) {

  double mu = 0; double delta = 0; std::vector<double> st(4);

  delta = shape / sqrt(1. + pow(shape, 2.));
  mu = delta * sqrt(df / M_PI) * exp(lgamma((df - 1.) / 2.) - lgamma(df / 2.));

  st[0] = mu;
  st[1] = df / (df - 2.) - pow(mu, 2.); // * escala ^ 2
  st[2] = mu * (df * (3. - pow(delta, 2.)) / (df - 3.) - 3. * df / (df - 2.) + 2. * pow(mu, 2.));
  st[3] = (3. * pow(df, 2.) / ((df - 2) * (df - 4.)) -
    4. * pow(mu, 2.) * df * (3. - pow(delta, 2.)) / (df - 3.) +
    6. * pow(mu, 2.) * df / (df - 2.) - 3. * pow(mu, 4.)) - 3. * pow(st[1], 2.);

  // FALTA ESTA PARTE COMENTADA PARA QUE LA FUNCION ESTE IGUAL QUE EN LA LIBRERIA DE R 'sn'
  // ASI ESTA IGUAL QUE EN EL EXCEL, ESTO NO AFECTA YA QUE SE IMPONE QUE EN LA OPTIMIZACION
  // QUE LOCATION = 0 y ESCALA = 1

  // st[0] = st[0] * escala + location;
  // st[1] = st[1] * pow(escala, 2);
  // st[2] = st[2] * pow(escala, 3);
  // st[3] = st[3] * pow(escala, 4);

  return st;
}


// Funcion de minimizacion para el ajuste de la skewt
// [[Rcpp::export]]
double minfunc(std::vector<double> par, std::vector<double> x, int n_days) {

  // if (par.size() != 2) stop("par vector must have two elements");
  int l = x.size();
  double k1o = k1(x);

  std::vector<double> C_X(l);

  for (int i = 0; i < l; i++) {

    C_X[i] = x[i] - k1o;

  }

  double k2c = k2(C_X) * n_days;
  double k3c = k3(C_X) * n_days;
  double k4c = k4(C_X) * n_days;


  double eg1 = k3c / pow(k2c, 1.5);
  double eg2 = k4c / pow(k2c, 2.);

  double location = 0;
  double escala = 1;
  double shape = par[0];
  double expDf = par[1];
  double df = exp(expDf) + 4.;

  std::vector<double> st(4);
  st = st_cumulants(location, escala, shape, df);

  double pg1 = st[2] / pow(st[1], 1.5);
  double pg2 = st[3] / pow(st[1], 2.);

  double min = 0;

  min = pow(std::abs(eg1 - pg1), 1.5) / (1 + std::abs(pg1)) + pow(std::abs(eg2 - pg2), 1.5) / (1 + pg2);

  return min;
}


// [[Rcpp::export]]
double k4m(std::vector<double> k, double s) {

  double k4m = 0;

  k4m = k[0] * s +
    (k[1] / 2) * std::pow(s, 2.) +
    (k[2] / 6) * std::pow(s, 3.) +
    (k[3] / 24) * std::pow(s, 4.) +
    (k[4] / 120) * std::pow(s, 5.) +
    (k[5] / 720) * std::pow(s, 6.);

  return k4m;

}


// [[Rcpp::export]]
double k4m_1st(std::vector<double> k, double s) {

  double k4m_1st = 0;

  k4m_1st = k[0] +
    k[1] * s +
    (k[2] / 2) * std::pow(s, 2) +
    (k[3] / 6) * std::pow(s, 3) +
    (k[4] / 24) * std::pow(s, 4) +
    (k[5] / 120) * std::pow(s, 5);

  return k4m_1st;

}


// [[Rcpp::export]]
double k4m_2nd(std::vector<double> k, double s) {

  double k4m_2nd = 0;

  k4m_2nd = k[1] +
    k[2] * s +
    (k[3] / 2) * std::pow(s, 2) +
    (k[4] / 6 )* std::pow(s, 3) +
    (k[5] / 24) * std::pow(s, 4);

  return k4m_2nd;

}
