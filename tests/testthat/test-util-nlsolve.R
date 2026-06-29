# util_nlsolve (issue #27): thin wrapper over nleqslv::nleqslv and BB::dfsane.
# Previously untested because it was only reached through the equilibrium_solve_*
# demography path. Tested here directly on a small system with a known root, plus
# the non-convergence error path. Pure numerics -- no model, no SCM.

# System: x^2 + y^2 = 4, x*y = 1. The root near (1.5, 1) is
# (x, y) = ((1+sqrt(3))/sqrt(2), (sqrt(3)-1)/sqrt(2)) ~= (1.93185, 0.51764).
sys <- function(z) c(z[1]^2 + z[2]^2 - 4, z[1] * z[2] - 1)
root <- c((1 + sqrt(3)) / sqrt(2), (sqrt(3) - 1) / sqrt(2))

test_that("util_nlsolve finds a known root with both solvers", {
  for (s in c("nleqslv", "dfsane")) {
    r <- util_nlsolve(c(1.5, 1), sys, solver = s)
    expect_equal(as.numeric(r), root, tolerance = 1e-4, info = s)
    expect_lt(max(abs(sys(r))), 1e-5)            # residual ~0
    expect_true(attr(r, "converged"), info = s)
    expect_equal(attr(r, "solver"), s)
  }
})

test_that("util_nlsolve respects the requested tolerance", {
  loose <- util_nlsolve(c(1.5, 1), sys, tol = 1e-3)
  tight <- util_nlsolve(c(1.5, 1), sys, tol = 1e-10)
  expect_lte(max(abs(sys(tight))), max(abs(sys(loose))) + 1e-12)
})

test_that("util_nlsolve rejects an unknown solver", {
  expect_error(util_nlsolve(c(1, 1), sys, solver = "newton"), "should be one of")
})

test_that("util_nlsolve errors when the system has no solution", {
  # x^2 = -1, y^2 = -1 has no real root; the solver cannot converge.
  no_root <- function(z) c(z[1]^2 + 1, z[2]^2 + 1)
  expect_error(util_nlsolve(c(0, 0), no_root, maxit = 5), "Solver has likely failed")
})
