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

// ---------------------------------------------------------------------------
// Dieckmann & Doebeli 1999 competition model
// ---------------------------------------------------------------------------
// "On the origin of species by sympatric speciation", Nature 400:354-357.
// Continuous-time logistic competition for a Gaussian resource:
//
//   K(x)   = K0 * exp(-(x - x0)^2 / (2 sigma_K^2))   carrying capacity
//   C(d)   = exp(-d^2 / (2 sigma_C^2))               competition kernel (C(0)=1)
//   single-resident equilibrium:  N* = K(x)
//   many-resident equilibrium:    solve A N = K, A_ij = C(x_i - x_j)
//   invasion fitness (a rate, not a ratio; ~0 at the resident):
//       s(y) = r * (1 - sum_i N_i C(y - x_i) / K(y))
//
// Analytic oracles: singular strategy x* = x0; it is a branching point
// (fitness minimum) iff sigma_C < sigma_K, an ESS (maximum) iff sigma_C > sigma_K.

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
//' Single resident: N* = K(x). Many residents: solve the linear system
//' A N = K with A_ij = C(x_i - x_j) (negative densities, i.e. strategies that
//' cannot coexist, are clamped to zero).
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
  std::vector<std::vector<double>> A(nr, std::vector<double>(nr + 1));
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

// ---------------------------------------------------------------------------
// Geritz, Kisdi, Meszena & Metz 1998 soft-selection model
// ---------------------------------------------------------------------------
// "Evolutionarily singular strategies and the adaptive growth and branching of
// the evolutionary tree", Evol. Ecol. 12:35-57 (the worked example: a haploid
// Levene soft-selection model). m patches, each with capacity K_j and a Gaussian
// within-patch survival peaked at optimum mu_j:
//
//   f_j(x) = exp(-(x - mu_j)^2 / (2 sigma^2))         (amplitude cancels)
//   resident equilibrium (Eq 19):  N'_i = sum_j K_j f_j(x_i) N_i / (sum_h f_j(x_h) N_h)
//   invasion fitness (Eq 21):  S(y) = log sum_j K_j f_j(y) / (sum_h f_j(x_h) N_h)
//   (for a single resident this reduces to Eq B1:
//        S(y) = log sum_j c_j f_j(y)/f_j(x),  c_j = K_j / sum K)
//
// Analytic oracles (symmetric m=3, mu=(-d,0,d), equal K): singular strategy
// x* = sum c_j mu_j = 0; convergence stable for all d/sigma; a branching point
// (lacks ESS) iff d/sigma > sqrt(3/2) ~= 1.2247, a CSS below that.

namespace geritz {

static inline double f_of(double x, double mu, double sigma) {
  double d = x - mu;
  return std::exp(-(d * d) / (2.0 * sigma * sigma));
}

} // namespace geritz

//' Geritz 1998 soft-selection model: log invasion fitness of mutants
//'
//' @param x_mut numeric vector of mutant trait values
//' @param x_res numeric vector of resident trait values
//' @param n_res numeric vector of resident equilibrium abundances
//' @param pars list with sigma, mu (patch optima), K (patch capacities)
//' @return numeric vector of log invasion fitness (=0 for resident at equilibrium)
//' @keywords internal
// [[Rcpp::export]]
NumericVector geritz_log_fitness(NumericVector x_mut, NumericVector x_res,
                                 NumericVector n_res, List pars) {
  double sigma = pars["sigma"];
  NumericVector mu = pars["mu"], Kp = pars["K"];
  int m = mu.size();
  int nm = x_mut.size();
  int nr = x_res.size();
  NumericVector out(nm);

  if (nr == 0) {
    // fundamental fitness of a lone strategy: log sum_j c_j f_j(y)
    double Ktot = 0.0;
    for (int j = 0; j < m; j++) Ktot += Kp[j];
    for (int i = 0; i < nm; i++) {
      double s = 0.0;
      for (int j = 0; j < m; j++)
        s += (Kp[j] / Ktot) * geritz::f_of(x_mut[i], mu[j], sigma);
      out[i] = std::log(s);
    }
    return out;
  }

  // patch denominators D_j = sum_h f_j(x_h) N_h
  std::vector<double> D(m, 0.0);
  for (int j = 0; j < m; j++)
    for (int h = 0; h < nr; h++)
      D[j] += geritz::f_of(x_res[h], mu[j], sigma) * n_res[h];

  for (int i = 0; i < nm; i++) {
    double s = 0.0;
    for (int j = 0; j < m; j++)
      s += Kp[j] * geritz::f_of(x_mut[i], mu[j], sigma) / D[j];
    out[i] = std::log(s);
  }
  return out;
}

//' Geritz 1998 soft-selection model: resident equilibrium abundances
//'
//' Iterates the soft-selection recursion (Eq 19) to its fixed point; total
//' abundance is conserved at sum(K). Single resident: N* = sum(K).
//'
//' @param x_res numeric vector of resident trait values
//' @param pars list with sigma, mu (patch optima), K (patch capacities)
//' @param max_iter maximum fixed-point iterations
//' @param eps convergence tolerance
//' @return numeric vector of equilibrium abundances
//' @keywords internal
// [[Rcpp::export]]
NumericVector geritz_equilibrium(NumericVector x_res, List pars,
                                 int max_iter = 5000, double eps = 1e-12) {
  double sigma = pars["sigma"];
  NumericVector mu = pars["mu"], Kp = pars["K"];
  int m = mu.size();
  int nr = x_res.size();
  NumericVector n(nr);
  if (nr == 0) return n;
  double Ktot = 0.0;
  for (int j = 0; j < m; j++) Ktot += Kp[j];
  if (nr == 1) {
    n[0] = Ktot;
    return n;
  }
  for (int h = 0; h < nr; h++) n[h] = Ktot / nr; // initial guess

  // precompute f_j(x_h)
  std::vector<std::vector<double>> F(m, std::vector<double>(nr));
  for (int j = 0; j < m; j++)
    for (int h = 0; h < nr; h++) F[j][h] = geritz::f_of(x_res[h], mu[j], sigma);

  for (int it = 0; it < max_iter; it++) {
    std::vector<double> D(m, 0.0);
    for (int j = 0; j < m; j++)
      for (int h = 0; h < nr; h++) D[j] += F[j][h] * n[h];
    double maxchange = 0.0;
    NumericVector nn(nr);
    for (int h = 0; h < nr; h++) {
      double v = 0.0;
      for (int j = 0; j < m; j++) v += Kp[j] * F[j][h] * n[h] / D[j];
      nn[h] = v;
      maxchange = std::max(maxchange, std::abs(nn[h] - n[h]));
    }
    n = nn;
    if (maxchange < eps) break;
  }
  return n;
}
