# Tests for community_solve_singularity_1D — the 1D evolutionary attractor
# (selection gradient == 0). Integration tests: each runs the SCM.

test_that("community_solve_singularity_1D finds the 1D attractor", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_solve_singularity_1D(bounds = c(0.13, 0.16), tol = 1e-3)

  expect_s3_class(comm, "community")
  expect_length(comm, 1L)                       # single resident at the singularity

  # the attractor sits inside the search bracket (reference ~0.1417 for current
  # plant 2.0.0.9001)
  x <- as.numeric(comm$traits[, "lma"])
  expect_gt(x, 0.13)
  expect_lt(x, 0.16)
  expect_equal(x, 0.1417, tolerance = 1e-2)

  # the selection gradient vanishes at the singularity. The root tolerance is on
  # the trait, so the gradient is small relative to its scale (|grad| reaches
  # hundreds away from the attractor) rather than exactly zero.
  expect_true(is.finite(comm$selection_gradient))
  expect_lt(abs(comm$selection_gradient), 50)
})

test_that("community_solve_singularity_1D handles non-bracketing bounds", {
  comm0 <- community_start(bounds(lma = c(0.01, 2)),
                           model_support = assembly_model_support())

  # bounds entirely above the attractor (~0.1417): the gradient is negative at
  # both ends, so there is no root to bracket. With edge_ok = TRUE (default)
  # this warns and returns the lower edge; with edge_ok = FALSE it errors.
  expect_warning(
    edge <- community_solve_singularity_1D(comm0, bounds = c(0.2, 0.5),
                                           tol = 1e-3),
    "Bounds do not include attractor")
  expect_equal(as.numeric(edge$traits[, "lma"]), 0.2)

  expect_error(
    community_solve_singularity_1D(comm0, bounds = c(0.2, 0.5), tol = 1e-3,
                                   edge_ok = FALSE),
    "Bounds do not include attractor")
})
