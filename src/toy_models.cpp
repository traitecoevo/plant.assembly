// Fast toy adaptive-dynamics models for developing and validating regnans's
// assembly / attractor machinery against models with known analytic answers
// (issue #33). Each model exposes an invasion-fitness function and a
// resident-equilibrium solver. Kept deliberately simple and dependency-free.
//
// Models implemented here:
//   - bird : migratory-bird arrival time (Johansson & Jonzen 2012;
//            Brannstrom, Johansson & von Festenberg 2013, Games 4:304-328, sec 4)
//
// Convention shared by all models: invasion-fitness functions return the
// LOG invasion fitness (so the resident's own value is ~0, matching the
// plant backend which returns log net-reproduction ratios).

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// ---------------------------------------------------------------------------
// Bird arrival-time model
// ---------------------------------------------------------------------------
// Trait x = arrival time. Birds compete for K territories; early arrival raises
// competitive ability C(x) = exp(-a x), reproduction R(x) is Gaussian about the
// seasonal optimum x_opt, survival p in (0,1).
//
//   R(x) = R0 * exp(-(x - x_opt)^2 / (2 sigma^2))
//   C(x) = exp(-a x)
//   resident demography (discrete):  n_{t+1} = K R(x) + p n_t
//   single-resident equilibrium:     n* = K R(x) / (1 - p)
//   invasion fitness of mutant x':    w = K R(x') C(x') / (sum_j n_j C(x_j)) + p
//   (returns log w; for one resident at n* this reduces to
//    w = (1-p) R(x')C(x') / (R(x)C(x)) + p, and w(x,x) = 1 exactly.)
//
// Analytic oracle: singular strategy x* = x_opt - a*sigma^2 (a CSS).

namespace bird {

static inline double R_of(double x, double R0, double x_opt, double sigma) {
  double d = x - x_opt;
  return R0 * std::exp(-(d * d) / (2.0 * sigma * sigma));
}

static inline double C_of(double x, double a) {
  return std::exp(-a * x);
}

} // namespace bird

//' Bird model: log invasion fitness of mutants
//'
//' @param x_mut numeric vector of mutant trait values (arrival times)
//' @param x_res numeric vector of resident trait values
//' @param n_res numeric vector of resident equilibrium densities
//' @param pars list with a, x_opt, sigma, R0, K, p
//' @return numeric vector of log invasion fitness, one per mutant
//' @keywords internal
// [[Rcpp::export]]
NumericVector bird_log_fitness(NumericVector x_mut, NumericVector x_res,
                               NumericVector n_res, List pars) {
  double a = pars["a"], x_opt = pars["x_opt"], sigma = pars["sigma"],
         R0 = pars["R0"], K = pars["K"], p = pars["p"];
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);

  if (nr == 0) {
    // Fundamental fitness of a lone strategy: log of its equilibrium
    // reproduction. Positive everywhere R > 0 (bird model is always viable),
    // so this is only a sensible-default fallback, not used on the main path.
    for (int i = 0; i < nm; i++) {
      out[i] = std::log(K * bird::R_of(x_mut[i], R0, x_opt, sigma) / (1.0 - p));
    }
    return out;
  }

  double denom = 0.0;
  for (int j = 0; j < nr; j++) denom += n_res[j] * bird::C_of(x_res[j], a);

  for (int i = 0; i < nm; i++) {
    double w = K * bird::R_of(x_mut[i], R0, x_opt, sigma) *
                   bird::C_of(x_mut[i], a) / denom + p;
    out[i] = std::log(w);
  }
  return out;
}

//' Bird model: resident demographic equilibrium densities
//'
//' Single resident: closed form n* = K R(x)/(1-p). Multiple residents: iterate
//' the territory-competition recursion to its fixed point (the bird model has a
//' single limiting resource, so this generically resolves to competitive
//' exclusion of all but the strategy maximising R(x)C(x)).
//'
//' @param x_res numeric vector of resident trait values
//' @param pars list with a, x_opt, sigma, R0, K, p
//' @param max_iter maximum fixed-point iterations (multi-resident case)
//' @param eps convergence tolerance (multi-resident case)
//' @return numeric vector of equilibrium densities, one per resident
//' @keywords internal
// [[Rcpp::export]]
NumericVector bird_equilibrium(NumericVector x_res, List pars,
                               int max_iter = 5000, double eps = 1e-12) {
  double a = pars["a"], x_opt = pars["x_opt"], sigma = pars["sigma"],
         R0 = pars["R0"], K = pars["K"], p = pars["p"];
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) {
    n[0] = K * bird::R_of(x_res[0], R0, x_opt, sigma) / (1.0 - p);
    return n;
  }

  std::vector<double> Rv(nr), Cv(nr);
  for (int j = 0; j < nr; j++) {
    Rv[j] = bird::R_of(x_res[j], R0, x_opt, sigma);
    Cv[j] = bird::C_of(x_res[j], a);
    n[j] = K * Rv[j] / (1.0 - p) / nr; // initial guess
  }

  for (int it = 0; it < max_iter; it++) {
    double denom = 0.0;
    for (int j = 0; j < nr; j++) denom += n[j] * Cv[j];
    double maxchange = 0.0;
    NumericVector nn(nr);
    for (int j = 0; j < nr; j++) {
      nn[j] = K * (n[j] * Cv[j] / denom) * Rv[j] + p * n[j];
      maxchange = std::max(maxchange, std::abs(nn[j] - n[j]));
    }
    n = nn;
    if (maxchange < eps) break;
  }
  return n;
}
