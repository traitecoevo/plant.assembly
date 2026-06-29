# community_fitness_landscape, grid method (issue #27). The landscape machinery
# is model-agnostic, so it is exercised on the fast DD99 harness rather than the
# SCM. (The bayesopt/surrogate method is not covered here -- it pulls in
# mlr3mbo and is exercised separately.)

dd99_resident <- function(...) {
  community_start(bounds(x = c(-2, 2)),
                  harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
                  trait_scale = "linear",
                  fitness_control = list(method = "grid", ...))
}

test_that("community_fitness_landscape (grid) samples fitness over the bounds", {
  comm <- dd99_resident(n_evals = 12) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography() |>
    community_fitness_landscape()

  pts <- comm$fitness_points
  expect_s3_class(pts, "tbl_df")
  # first column is named by the trait, plus fitness + resident flag
  expect_equal(names(pts), c("x", "fitness", "resident"))

  # grid spans the bounds (augmented by the resident +/- offset points)
  expect_gte(nrow(pts), 12L)
  expect_gte(min(pts$x), -2)
  expect_lte(max(pts$x), 2)
  expect_true(all(is.finite(pts$fitness)))

  # the resident is flagged and (at equilibrium) has fitness ~0
  expect_true(any(pts$resident))
  expect_equal(pts$x[pts$resident], 0.5, tolerance = 1e-8)
  expect_lt(abs(pts$fitness[pts$resident]), 1e-6)
})

test_that("community_fitness_landscape solves demography if needed", {
  # No community_demography() call: the landscape function should solve it.
  comm <- community_fitness_landscape(dd99_resident(n_evals = 8))
  expect_true(is.function(comm$fitness_function))
  expect_s3_class(comm$fitness_points, "tbl_df")
})

test_that("community_fitness_landscape rejects an unknown method", {
  comm <- dd99_resident(n_evals = 8) |>
    community_add(trait_matrix(0.5, "x"), birth_rate = 100) |>
    community_demography()
  expect_error(community_fitness_landscape(comm, method = "nope"),
               "Unknown fitness landscape method")
})
