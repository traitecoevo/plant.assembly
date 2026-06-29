# Demographic equilibrium solvers (issue #27).
#
# community_demography() dispatches on demography_control$equilibrium_solver_name
# to one of five backends. Previously only the default (equilibrium_iteration)
# was tested, and only via the plant SCM (seconds per solve); the
# equilibrium_solve_* paths -- which are the only callers of util_nlsolve -- were
# untested against the community object (see #27 "blocked/deferred").
#
# Here we drive all five through the DD99 toy harness, whose multi-resident
# equilibrium has an analytic answer (the competition linear-solve), so every
# solver can be checked against a known target in milliseconds.
#
# NOTE: the explicit (toy) harness returns its equilibrium directly from the
# demography runner, so for DD99 the solvers do not have to *iterate* to a fixed
# point the way the plant SCM does. What these tests verify is the dispatch, the
# util_nlsolve integration, and the cleanup write-back -- all against the
# analytic oracle. The genuinely-iterative convergence of the default solver is
# covered on the real SCM in test-plant-smoke.R.

dd99_pars <- list(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4)

solve_with <- function(solver, x = c(-0.5, 0.5)) {
  community_start(
    bounds(x = c(-2, 2)),
    harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
    demography_control = demographic_step_control(
      list(equilibrium_solver_name = solver))) |>
    community_add(trait_matrix(x, "x"), birth_rate = rep(100, length(x))) |>
    community_demography()
}

test_that("every equilibrium solver recovers the analytic DD99 equilibrium", {
  target <- dd99_equilibrium(c(-0.5, 0.5), dd99_pars)
  solvers <- c("equilibrium_iteration", "single_step", "equilibrium_hybrid",
               "equilibrium_solve_nleqslv", "equilibrium_solve_dfsane")
  for (s in solvers) {
    comm <- solve_with(s)
    expect_true(attr(comm, "converged"), info = s)
    expect_equal(as.numeric(comm$birth_rate), target, tolerance = 1e-5, info = s)
    # at the equilibrium each resident's invasion fitness is ~0
    expect_equal(comm$resident_fitness, c(0, 0), tolerance = 1e-6, info = s)
  }
})

test_that("the nleqslv and dfsane solvers agree with the default iteration", {
  ref <- as.numeric(solve_with("equilibrium_iteration")$birth_rate)
  expect_equal(as.numeric(solve_with("equilibrium_solve_nleqslv")$birth_rate),
               ref, tolerance = 1e-5)
  expect_equal(as.numeric(solve_with("equilibrium_solve_dfsane")$birth_rate),
               ref, tolerance = 1e-5)
})

test_that("community_demography rejects an unknown solver", {
  comm <- community_start(
    bounds(x = c(-2, 2)), harness = harness_dd99(),
    demography_control = demographic_step_control(
      list(equilibrium_solver_name = "not_a_solver"))) |>
    community_add(trait_matrix(0, "x"), birth_rate = 100)
  expect_error(community_demography(comm), "Unknown solver")
})

test_that("a three-resident community also solves to the analytic equilibrium", {
  x <- c(-1, 0, 1)
  target <- dd99_equilibrium(x, dd99_pars)
  comm <- solve_with("equilibrium_solve_nleqslv", x = x)
  expect_equal(as.numeric(comm$birth_rate), target, tolerance = 1e-5)
  expect_equal(comm$resident_fitness, rep(0, 3), tolerance = 1e-6)
})
