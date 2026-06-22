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

# Viable (fundamental) fitness, reimplemented on the community machinery.
# Replaces the old test of plant's removed fundamental_fitness()/viable_fitness()
# /fitness_landscape()/max_fitness(); those low-level plant functions no longer
# exist (plant #388). The fundamental niche is now computed via an empty
# community's fitness_function.
test_that("community_viable_fitness_1D finds the positive-fitness region", {
  ms <- list(p = plant_default_assembly_pars(max_patch_lifetime = 30),
             plant_control = plant_default_assembly_control())
  comm <- community_start(bounds(lma = c(0.01, 2)), model_support = ms)

  vb <- community_viable_fitness_1D(comm)
  expect_true(is.matrix(vb))
  expect_equal(rownames(vb), "lma")
  expect_equal(colnames(vb), c("lower", "upper"))

  # the viable interval sits strictly inside the search bounds
  expect_gt(vb[, "lower"], 0.01)
  expect_lt(vb[, "upper"], 2)
  expect_lt(vb[, "lower"], vb[, "upper"])

  # reference values for current plant (2.0.0.9001)
  expect_equal(unname(vb[, "lower"]), 0.01309, tolerance = 1e-3)
  expect_equal(unname(vb[, "upper"]), 1.64417, tolerance = 1e-3)

  # fitness is positive inside the interval and ~0 at the edges
  fitness <- plant_community_update_fitness_function(comm)$fitness_function
  mid <- exp(mean(log(vb)))
  expect_gt(fitness(mid), 1)                    # clearly positive inside
  expect_lt(abs(fitness(vb[, "lower"])), 0.01)  # ~0 at the edges (root tol is
  expect_lt(abs(fitness(vb[, "upper"])), 0.01)  # on the trait, not on fitness)
})

# max_fitness / max_growth_rate, reimplemented on the community machinery
# (replacing plant's removed plant-level versions; resolves the previous
# duplicate max_fitness definitions).
test_that("max_growth_rate evaluates community fitness at trait values", {
  ms <- list(p = plant_default_assembly_pars(max_patch_lifetime = 30),
             plant_control = plant_default_assembly_control())
  comm <- community_start(bounds(lma = c(0.01, 2)), model_support = ms)

  g <- max_growth_rate(comm, c(0.05, 0.0825, 0.2))
  expect_length(g, 3L)
  expect_true(all(is.finite(g)))
  # agrees with the community's fitness_function directly
  ff <- plant_community_update_fitness_function(comm)$fitness_function
  expect_equal(g, ff(c(0.05, 0.0825, 0.2)), tolerance = 1e-8)
})

test_that("max_fitness finds the fitness peak within bounds", {
  ms <- list(p = plant_default_assembly_pars(max_patch_lifetime = 30),
             plant_control = plant_default_assembly_control())
  comm <- community_start(bounds(lma = c(0.01, 2)), model_support = ms)

  mx <- max_fitness(comm)
  expect_equal(names(mx), "lma")
  w <- attr(mx, "fitness")
  expect_true(is.finite(w))

  # the maximiser sits inside the bounds and its fitness is the largest seen
  expect_gt(as.numeric(mx), 0.01)
  expect_lt(as.numeric(mx), 2)
  grid <- max_growth_rate(comm, seq_log_range(c(0.01, 2), 15))
  expect_gte(w, max(grid) - 1e-3)
})
