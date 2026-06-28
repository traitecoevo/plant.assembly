// GM99 -- Geritz, van der Meijden & Metz 1999 seed-size safe-site model
// "Evolutionary Dynamics of Seed Size and Seedling Competitive Ability."
//
// Plants compete for safe sites; seeds arrive by a Poisson process; seedlings
// undergo size-asymmetric lottery competition within each site.
//   trait x = seed size, 0 < x <= R (R = resources per site)
//   pre-competitive survival   s(x) = max(0, 1 - 2 exp(-beta x))
//   competitive ability        c(x) = exp(alpha x)
//   establishment of mutant x' against residents {x_j} at densities {N_j}:
//     g = E[ c(x') / (c(x') + sum_j k_j c(x_j)) ],  k_j ~ Poisson(N_j) indep.
//   invasion fitness (lifetime reproductive ratio, = 1 at equilibrium):
//     W = (R/x') s(x') g          (gm99_fitness returns log W; resident ~0)
//   resident equilibrium: W_i = 1 for every resident; single resident solves
//     (R/x) s(x) (1 - exp(-N))/N = 1   (exists iff (R/x) s(x) > 1).
//
// Only the products alpha*R and beta*R matter. No closed-form singular strategy
// (numerical); branching diagnosed by the sign of the invasion-fitness second
// derivative. This is the model implemented in Daniel's MATLAB.
// Distinct from the GK98 1998 soft-selection model.

#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace gm99 {

static inline double s_of(double x, double beta) {
  double v = 1.0 - 2.0 * std::exp(-beta * x);
  return v > 0.0 ? v : 0.0;
}
static inline double c_of(double x, double alpha) {
  return std::exp(alpha * x);
}

// E[ c_mut / (c_mut + sum_j k_j cres_j) ] over independent k_j ~ Poisson(N_j),
// using a numerically stable Poisson pmf recurrence (p_0 = e^-N, p_{k+1} =
// p_k N/(k+1)) and an odometer over the per-resident truncations.
static double estab(double c_mut, const std::vector<double>& cres,
                    const std::vector<double>& N) {
  int n = cres.size();
  if (n == 0) return 1.0;
  std::vector<std::vector<double> > P(n);
  std::vector<int> Kmax(n);
  for (int j = 0; j < n; j++) {
    int kmax = (int)std::ceil(N[j] + 10.0 * std::sqrt(N[j] + 1.0) + 20.0);
    P[j].resize(kmax + 1);
    double pk = std::exp(-N[j]);
    for (int k = 0; k <= kmax; k++) { P[j][k] = pk; pk *= N[j] / (k + 1); }
    Kmax[j] = kmax;
  }
  std::vector<int> k(n, 0);
  double sum = 0.0;
  while (true) {
    double w = 1.0, denom = c_mut;
    for (int j = 0; j < n; j++) { w *= P[j][k[j]]; denom += k[j] * cres[j]; }
    sum += w * c_mut / denom;
    int j = 0;
    while (j < n) { if (++k[j] <= Kmax[j]) break; k[j] = 0; j++; }
    if (j == n) break;
  }
  return sum;
}

} // namespace gm99

//' GM99 seed-size model: log invasion fitness of mutants
//'
//' @param x_mut numeric vector of mutant seed sizes
//' @param x_res numeric vector of resident seed sizes
//' @param n_res numeric vector of resident equilibrium densities (seeds/site)
//' @param pars list with R, alpha, beta
//' @return numeric vector of log invasion fitness (=0 for resident at equilibrium;
//'   -Inf for non-viable mutants where s(x) = 0)
//' @keywords internal
// [[Rcpp::export]]
NumericVector gm99_fitness(NumericVector x_mut, NumericVector x_res,
                           NumericVector n_res, List pars) {
  double R = pars["R"], alpha = pars["alpha"], beta = pars["beta"];
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);
  std::vector<double> cres(nr), N(nr);
  for (int j = 0; j < nr; j++) { cres[j] = gm99::c_of(x_res[j], alpha); N[j] = n_res[j]; }
  for (int i = 0; i < nm; i++) {
    double mp = x_mut[i];
    double s = gm99::s_of(mp, beta);
    if (s <= 0.0 || mp <= 0.0) { out[i] = R_NegInf; continue; }
    double seeds = R / mp * s;
    double g = (nr == 0) ? 1.0 : gm99::estab(gm99::c_of(mp, alpha), cres, N);
    out[i] = std::log(seeds * g);
  }
  return out;
}

//' GM99 seed-size model: resident equilibrium densities (seeds/site)
//'
//' Single resident: solve (R/x) s(x) (1 - exp(-N))/N = 1 by bisection (N* = 0 if
//' the strategy is non-viable). Many residents: iterate the population recursion
//' N_i <- N_i W_i to its fixed point.
//'
//' @param x_res numeric vector of resident seed sizes
//' @param pars list with R, alpha, beta
//' @param max_iter maximum iterations (multi-resident case)
//' @param eps convergence tolerance (multi-resident case)
//' @return numeric vector of equilibrium densities
//' @keywords internal
// [[Rcpp::export]]
NumericVector gm99_equilibrium(NumericVector x_res, List pars,
                               int max_iter = 10000, double eps = 1e-12) {
  double R = pars["R"], alpha = pars["alpha"], beta = pars["beta"];
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  if (nr == 1) {
    double m = x_res[0], s = gm99::s_of(m, beta);
    double seeds = (s > 0.0 && m > 0.0) ? R / m * s : 0.0;
    if (seeds <= 1.0) { n[0] = 0.0; return n; }
    // (R/x) s(x) (1 - e^-N)/N = 1 is monotone decreasing in N; bisect [0, seeds]
    double lo = 1e-12, hi = seeds;
    for (int it = 0; it < 100; it++) {
      double mid = 0.5 * (lo + hi);
      double f = seeds * (1.0 - std::exp(-mid)) / mid - 1.0;
      if (f > 0.0) lo = mid; else hi = mid;
    }
    n[0] = 0.5 * (lo + hi);
    return n;
  }
  // multi-resident: iterate N_i <- N_i * W_i
  std::vector<double> cres(nr);
  for (int j = 0; j < nr; j++) {
    cres[j] = gm99::c_of(x_res[j], alpha);
    NumericVector xj(1); xj[0] = x_res[j];
    n[j] = gm99_equilibrium(xj, pars)[0];
    if (n[j] <= 0.0) n[j] = 0.01;
  }
  for (int it = 0; it < max_iter; it++) {
    std::vector<double> Nv(nr);
    for (int j = 0; j < nr; j++) Nv[j] = n[j];
    double maxchange = 0.0;
    NumericVector nn(nr);
    for (int i = 0; i < nr; i++) {
      double m = x_res[i], s = gm99::s_of(m, beta);
      double seeds = (s > 0.0 && m > 0.0) ? R / m * s : 0.0;
      double W = seeds * gm99::estab(cres[i], cres, Nv);
      double newN = n[i] * W;
      if (newN < 0.0) newN = 0.0;
      maxchange = std::max(maxchange, std::abs(newN - n[i]));
      nn[i] = newN;
    }
    n = nn;
    if (maxchange < eps) break;
  }
  return n;
}
