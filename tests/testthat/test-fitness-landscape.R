# Tests for community_fitness_landscape (grid method). The bayesopt/surrogate
# method is not covered here (it pulls in mlr3mbo and is exercised separately).

test_that("community_fitness_landscape (grid) samples fitness over the bounds", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          fitness_control = list(method = "grid", n_evals = 12),
                          model_support = assembly_model_support()) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography() |>
    community_fitness_landscape()

  pts <- comm$fitness_points
  expect_s3_class(pts, "tbl_df")
  # first column is named by the trait, plus fitness + resident flag
  expect_equal(names(pts), c("lma", "fitness", "resident"))

  # grid spans the bounds, augmented by the resident +/- offset points
  expect_gte(nrow(pts), 12L)
  expect_gte(min(pts$lma), 0.01 * 0.99)
  expect_lte(max(pts$lma), 2)
  expect_true(all(is.finite(pts$fitness)))

  # the resident is flagged and (at equilibrium) has fitness ~0
  expect_true(any(pts$resident))
  expect_equal(pts$lma[pts$resident], 0.0825, tolerance = 1e-8)
  expect_lt(abs(pts$fitness[pts$resident]), 1e-3)
})

test_that("community_fitness_landscape solves demography if needed", {
  # No community_demography() call: the landscape function should solve it.
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          fitness_control = list(method = "grid", n_evals = 8),
                          model_support = assembly_model_support())
  comm <- community_fitness_landscape(comm)
  expect_true(is.function(comm$fitness_function))
  expect_s3_class(comm$fitness_points, "tbl_df")
})

test_that("community_fitness_landscape rejects an unknown method", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          fitness_control = list(method = "grid", n_evals = 8),
                          model_support = assembly_model_support()) |>
    community_demography()
  expect_error(community_fitness_landscape(comm, method = "nope"),
               "Unknown fitness landscape method")
})
