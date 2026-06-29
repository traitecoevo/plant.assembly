# community_solve_singularity_1D edge behaviour (issue #27).
#
# The happy path (recovering the singular strategy) is covered analytically on
# DD99 in test-harness-dd99.R ("dd99 singular strategy is x0") and on the real
# SCM in test-plant-smoke.R. Here we cover the non-bracketing branches, which
# are model-agnostic and so run on the fast DD99 harness: DD99's singular point
# is x0 = 0, so a search bracket sitting entirely above it has a same-signed
# (negative) gradient at both ends and cannot bracket a root.

dd99_empty <- function() {
  community_start(bounds(x = c(-2, 2)),
                  harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4))
}

test_that("community_solve_singularity_1D warns and returns the edge when bounds miss the attractor", {
  expect_warning(
    edge <- community_solve_singularity_1D(dd99_empty(), bounds = c(0.5, 1.5),
                                           tol = 1e-3),
    "Bounds do not include attractor")
  expect_equal(as.numeric(edge$traits[, "x"]), 0.5)   # returns the lower edge
})

test_that("community_solve_singularity_1D errors on non-bracketing bounds when edge_ok = FALSE", {
  expect_error(
    community_solve_singularity_1D(dd99_empty(), bounds = c(0.5, 1.5),
                                   tol = 1e-3, edge_ok = FALSE),
    "Bounds do not include attractor")
})
