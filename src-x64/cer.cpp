#include <Rcpp.h>
#include <iostream>
using namespace Rcpp;


// [[Rcpp::export]]
std::vector<double> RcOut(double A, double B, double C, double D,
                          std::vector<double> PD, std::vector<double> LGD,
                          std::vector<double> weight, std::vector<double> EAD,
                          std::vector<double> CORR, std::vector<double> UL) {

  int l = PD.size();
  double sum = 0;
  double avg_pd = 0;
  double sum_ead = 0;
  double avg_ead = 0;
  double avg_cor = 0;
  std::vector<double> out(l);

  for (int i = 0; i < l; i++) {
    sum += weight[i];
  }

  for (int i = 0; i < l; i++) {
    avg_pd += PD[i] * weight[i];
  }

  avg_pd = avg_pd / sum;

  for (int i = 0; i < l; i++) {
    sum_ead += EAD[i] * weight[i];
  }

  avg_ead = sum_ead / sum;

  for (int i = 0; i < l; i++) {
    avg_cor += CORR[i] * weight[i];
  }

  avg_cor = pow(avg_cor / sum, 2);

  for (int i = 0; i < l; i++) {
    out[i] = UL[i] * (A + B * (EAD[i] / sum_ead) * log(EAD[i] / avg_ead) + C *
      log(PD[i] / avg_pd) + D * CORR[i] * CORR[i] * log(pow(CORR[i], 2) / avg_cor)) -
      PD[i] * LGD[i];
  }

  return out;

}


double median(std::vector<double> a) {

  int l = a.size();
  std::sort(&a[0], &a[l]);
  double median = l % 2 ? a[l / 2] : (a[l / 2 - 1] + a[l / 2]) / 2;
  return median;
}


// [[Rcpp::export]]
double resid(double A, double B, double C, double D,
             std::vector<double> PD, std::vector<double> LGD,
             std::vector<double> weight, std::vector<double> EAD,
             std::vector<double> CORR, std::vector<double> UL,
             std::vector<double> RcIn) {

  int l = weight.size();
  std::vector<double> out(l);
  double resid = 0;
  double sum = 0;

  for (int i = 0; i < l; i++) {
    sum += weight[i];
  }

  out = RcOut(A, B, C, D, PD, LGD, weight, EAD, CORR, UL);

  for (int i = 0; i < l; i++) {
    out[i] = pow(weight[i] * (out[i] - RcIn[i]), 2);
  }

  resid = sqrt(median(out) / sum);

  return resid;
}
