// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// mat
NumericMatrix mat(NumericVector FG, NumericVector FL, CharacterVector RU, CharacterVector col);
RcppExport SEXP _ECTools_mat(SEXP FGSEXP, SEXP FLSEXP, SEXP RUSEXP, SEXP colSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type FG(FGSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type FL(FLSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type RU(RUSEXP);
    Rcpp::traits::input_parameter< CharacterVector >::type col(colSEXP);
    rcpp_result_gen = Rcpp::wrap(mat(FG, FL, RU, col));
    return rcpp_result_gen;
END_RCPP
}
// fgyfl
List fgyfl(NumericVector x, double lim);
RcppExport SEXP _ECTools_fgyfl(SEXP xSEXP, SEXP limSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type lim(limSEXP);
    rcpp_result_gen = Rcpp::wrap(fgyfl(x, lim));
    return rcpp_result_gen;
END_RCPP
}
// k1
double k1(std::vector<double> x);
RcppExport SEXP _ECTools_k1(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(k1(x));
    return rcpp_result_gen;
END_RCPP
}
// k2
double k2(std::vector<double> x);
RcppExport SEXP _ECTools_k2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(k2(x));
    return rcpp_result_gen;
END_RCPP
}
// k3
double k3(std::vector<double> x);
RcppExport SEXP _ECTools_k3(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(k3(x));
    return rcpp_result_gen;
END_RCPP
}
// k4
double k4(std::vector<double> x);
RcppExport SEXP _ECTools_k4(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(k4(x));
    return rcpp_result_gen;
END_RCPP
}
// k5
double k5(std::vector<double> x);
RcppExport SEXP _ECTools_k5(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(k5(x));
    return rcpp_result_gen;
END_RCPP
}
// k6
double k6(std::vector<double> x);
RcppExport SEXP _ECTools_k6(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(k6(x));
    return rcpp_result_gen;
END_RCPP
}
// st_cumulants
std::vector<double> st_cumulants(double location, double escala, double shape, double df);
RcppExport SEXP _ECTools_st_cumulants(SEXP locationSEXP, SEXP escalaSEXP, SEXP shapeSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type location(locationSEXP);
    Rcpp::traits::input_parameter< double >::type escala(escalaSEXP);
    Rcpp::traits::input_parameter< double >::type shape(shapeSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(st_cumulants(location, escala, shape, df));
    return rcpp_result_gen;
END_RCPP
}
// minfunc
double minfunc(std::vector<double> par, std::vector<double> x, int n_days);
RcppExport SEXP _ECTools_minfunc(SEXP parSEXP, SEXP xSEXP, SEXP n_daysSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type par(parSEXP);
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type n_days(n_daysSEXP);
    rcpp_result_gen = Rcpp::wrap(minfunc(par, x, n_days));
    return rcpp_result_gen;
END_RCPP
}
// k4m
double k4m(std::vector<double> k, double s);
RcppExport SEXP _ECTools_k4m(SEXP kSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(k4m(k, s));
    return rcpp_result_gen;
END_RCPP
}
// k4m_1st
double k4m_1st(std::vector<double> k, double s);
RcppExport SEXP _ECTools_k4m_1st(SEXP kSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(k4m_1st(k, s));
    return rcpp_result_gen;
END_RCPP
}
// k4m_2nd
double k4m_2nd(std::vector<double> k, double s);
RcppExport SEXP _ECTools_k4m_2nd(SEXP kSEXP, SEXP sSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type k(kSEXP);
    Rcpp::traits::input_parameter< double >::type s(sSEXP);
    rcpp_result_gen = Rcpp::wrap(k4m_2nd(k, s));
    return rcpp_result_gen;
END_RCPP
}
