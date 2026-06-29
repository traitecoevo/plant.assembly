# community_check_for_inviable_strategies (issue #27): the model-agnostic
# connector that flags residents whose equilibrium density has collapsed below
# the extinction threshold, so the assembler can drop them. Exercised here on
# the DD99 explicit harness (the plant-specific path is covered in
# test-plant-smoke.R).

mk <- function(x, eps = NULL) {
  ctrl <- if (is.null(eps)) demographic_step_control()
          else demographic_step_control(list(equilibrium_extinct_birth_rate = eps))
  community_start(bounds(x = c(-2, 2)),
                  harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
                  demography_control = ctrl) |>
    community_add(trait_matrix(x, "x"), birth_rate = rep(100, length(x))) |>
    community_demography()
}

test_that("viable residents are not flagged for dropping", {
  op <- community_check_for_inviable_strategies(mk(c(0, 1.9)))
  expect_length(op, 2L)
  expect_true(all(op > 0))                       # positive equilibrium densities
  expect_false(any(attr(op, "drop")))            # none extinct
})

test_that("a resident below the extinction threshold is flagged", {
  # DD99 equilibrium density falls with distance from the optimum: at x = 1.9 it
  # is ~82, at x = 0 it is ~500. With the threshold raised to 100 only the
  # peripheral resident is flagged.
  op <- community_check_for_inviable_strategies(mk(c(0, 1.9), eps = 100))
  drop <- attr(op, "drop")
  expect_equal(drop, c(FALSE, TRUE))
  expect_lt(op[2], 100)
  expect_gt(op[1], 100)
})
