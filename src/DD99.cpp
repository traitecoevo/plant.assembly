// DD99 -- Dieckmann & Doebeli 1999 competition model
// "On the origin of species by sympatric speciation", Nature 400:354-357.
//
// Continuous-time logistic competition for a Gaussian resource:
//   K(x)   = K0 exp(-(x - x0)^2 / (2 sigma_K^2))    carrying capacity
//   C(d)   = exp(-d^2 / (2 sigma_C^2))              competition kernel (C(0)=1)
//   single-resident equilibrium:  N* = K(x)
//   many-resident equilibrium:    solve A N = K, A_ij = C(x_i - x_j)
//   invasion fitness (a per-capita growth RATE, not a ratio; ~0 at the resident):
//       s(y) = r (1 - sum_i N_i C(y - x_i) / K(y))
//
// Analytic oracles: singular strategy x* = x0; a branching point (fitness
// minimum) iff sigma_C < sigma_K, an ESS (maximum) iff sigma_C > sigma_K.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace dd99 {

static inline double K_of(double x, double K0, double x0, double sK) {
  double d = x - x0;
  return K0 * std::exp(-(d * d) / (2.0 * sK * sK));
}
static inline double C_of(double d, double sC) {
  return std::exp(-(d * d) / (2.0 * sC * sC));
}

} // namespace dd99

//' DD99 model: invasion fitness of mutants (a per-capita growth rate)
//'
//' @param x_mut numeric vector of mutant trait values
//' @param x_res numeric vector of resident trait values
//' @param n_res numeric vector of resident equilibrium densities
//' @param pars list with r, K0, x0, sigma_K, sigma_C
//' @return numeric vector of invasion fitness (=0 for a resident at equilibrium)
//' @keywords internal
// [[Rcpp::export]]
NumericVector dd99_fitness(NumericVector x_mut, NumericVector x_res,
                           NumericVector n_res, List pars) {
  double r = pars["r"], K0 = pars["K0"], x0 = pars["x0"],
         sK = pars["sigma_K"], sC = pars["sigma_C"];
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);
  for (int i = 0; i < nm; i++) {
    double comp = 0.0;
    for (int j = 0; j < nr; j++)
      comp += n_res[j] * dd99::C_of(x_mut[i] - x_res[j], sC);
    double Ky = dd99::K_of(x_mut[i], K0, x0, sK);
    out[i] = r * (1.0 - comp / Ky);
  }
  return out;
}

//' DD99 model: resident demographic equilibrium densities
//'
//' Single resident: N* = K(x). Many residents: solve the linear system A N = K
//' with A_ij = C(x_i - x_j) (negative densities, i.e. strategies that cannot
//' coexist, are clamped to zero).
//'
//' @param x_res numeric vector of resident trait values
//' @param pars list with r, K0, x0, sigma_K, sigma_C
//' @return numeric vector of equilibrium densities
//' @keywords internal
// [[Rcpp::export]]
NumericVector dd99_equilibrium(NumericVector x_res, List pars) {
  double K0 = pars["K0"], x0 = pars["x0"],
         sK = pars["sigma_K"], sC = pars["sigma_C"];
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) {
    n[0] = dd99::K_of(x_res[0], K0, x0, sK);
    return n;
  }
  // dense linear solve via Gaussian elimination with partial pivoting
  std::vector<std::vector<double> > A(nr, std::vector<double>(nr + 1));
  for (int i = 0; i < nr; i++) {
    for (int j = 0; j < nr; j++) A[i][j] = dd99::C_of(x_res[i] - x_res[j], sC);
    A[i][nr] = dd99::K_of(x_res[i], K0, x0, sK);
  }
  for (int col = 0; col < nr; col++) {
    int piv = col;
    for (int i = col + 1; i < nr; i++)
      if (std::abs(A[i][col]) > std::abs(A[piv][col])) piv = i;
    std::swap(A[col], A[piv]);
    double d = A[col][col];
    for (int j = col; j <= nr; j++) A[col][j] /= d;
    for (int i = 0; i < nr; i++) {
      if (i == col) continue;
      double f = A[i][col];
      for (int j = col; j <= nr; j++) A[i][j] -= f * A[col][j];
    }
  }
  for (int i = 0; i < nr; i++) n[i] = std::max(0.0, A[i][nr]);
  return n;
}
