test_that("positive_1d", {
  f <- function(x) -x^2 + 1
  tol <- 1e-8
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  tol <- 1e-3
  r <- positive_1d(f, 0, 0.1, tol=tol)
  expect_lt(r[[1]], -1 + tol)
  expect_gt(r[[2]], 1 - tol)

  expect_error(positive_1d(f, -2, 0.1, tol=tol), "no positive values")
})

test_that("bounds", {

  # First test object with infinite bounds
  expect_silent({
    bounds0 <- bounds_infinite("lma")
  })
  expect_true(is.matrix(bounds0))
  expect_equal(bounds0,
    matrix(c(-Inf, Inf), nrow=1, dimnames =list("lma",c("lower", "upper")))
  )

  # Manually created bounds

  expect_silent(
    bounds_1d <- bounds(lma=c(0.01, 10))
  )
  expect_silent(
    bounds_2d <- bounds(lma=c(0.01, 10), rho=c(1, 1000))
  )

  expect_true(is.matrix(bounds_1d))
  expect_equal(
    bounds_1d,
    matrix(c(0.01, 10), nrow = 1, dimnames = list("lma", c("lower", "upper")))
  )

  expect_true(is.matrix(bounds_2d))
  expect_equal(
    bounds_2d,
    matrix(c(0.01, 10, 1, 1000), byrow=TRUE, nrow = 2, dimnames = list(c("lma", "rho"), c("lower", "upper")))
  )

  expect_silent(
    check_bounds(bounds_1d)
  )
  expect_silent(
    check_bounds(bounds_2d)
  )

  expect_error(
    check_bounds(c(0,1))
  )
  expect_error(
    check_bounds(matrix(c(0.01, 10, 1, 1000)))
  )
  expect_error(
    check_bounds(bounds0, finite=TRUE)
  )
  expect_silent(
    check_bounds(bounds0)
  )


  # Points lie within bounds
  expect_silent(
    check_point(0.02, bounds_1d)
  )
  expect_error(
    check_point(0.001, bounds_1d)
  )
})

# NOTE (issue #27): the SCM-backed tests that used to live here --
# community_viable_fitness_1D, max_growth_rate, max_fitness -- are log-scale /
# fundamental-niche behaviours specific to the plant path, so they now live in
# test-plant-smoke.R (the consolidated, minimal set of genuine SCM tests). The
# model-agnostic pipeline is covered fast on DD99 elsewhere.
