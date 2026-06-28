# Community assembly with the toy harnesses (issue #33): going beyond pairwise
# invasibility plots to actually run the assembly algorithm, plus the linear
# trait-scale option and the multi-trait (2D) DD99 model.

test_that("trait_scale option is recorded and defaults to log", {
  expect_equal(community_start(bounds(x = c(-2, 2)),
                               harness = harness_dd99())$trait_scale, "log")
  expect_equal(community_start(bounds(x = c(-2, 2)), harness = harness_dd99(),
                               trait_scale = "linear")$trait_scale, "linear")
  tf <- community_trait_transform(list(trait_scale = "linear"))
  expect_equal(tf$fwd(-3), -3)             # identity on a negative value
  expect_equal(community_trait_transform(list(trait_scale = "log"))$fwd(exp(1)), 1)
})

test_that("DD99 assembles a packed community on a linear trait scale", {
  # traits centred at 0 -> needs the linear scale; sigma_C < sigma_K -> branching
  comm <- community_start(bounds(x = c(-3, 3)),
                          harness = harness_dd99(x0 = 0, sigma_K = 1, sigma_C = 0.4),
                          trait_scale = "linear",
                          fitness_control = list(method = "grid", n_evals = 100))
  a <- assembler_start(comm, assembler_control(list(birth_type = "maximum")))
  a <- assembler_run(a, 12)
  # limiting similarity: many coexisting species spread across trait space
  expect_gt(length(a$community), 4L)
  expect_true(all(a$community$birth_rate > 0))
  expect_true(min(a$community$traits) < -1 && max(a$community$traits) > 1)
})

test_that("GK98 branches to a dimorphism under assembly (seeded)", {
  comm <- community_start(bounds(x = c(-4, 4)),
                          harness = harness_gk98(d = 1.5, sigma = 1),
                          trait_scale = "linear",
                          fitness_control = list(method = "grid", n_evals = 100))
  a <- assembler_start(comm, assembler_control(list(birth_type = "maximum")))
  a <- assembler_set_traits(a, trait_matrix(0, "x"))  # soft selection: seed first type
  a <- assembler_run(a, 10)
  expect_gte(length(a$community), 2L)                  # branched
  # symmetric dimorphism: one type below 0, one above
  tr <- sort(as.numeric(a$community$traits))
  expect_lt(tr[1], 0); expect_gt(tr[length(tr)], 0)
})

test_that("DD99 multi-trait (2D) invasion fitness has its singular point at x0", {
  p <- list(r = 1, K0 = 500, x0 = c(0, 0), sigma_K = c(1, 1), sigma_C = c(0.5, 0.5))
  xr <- matrix(c(0, 0), nrow = 1)
  n <- dd99_nd_equilibrium(xr, p)
  expect_equal(as.numeric(n), 500, tolerance = 1e-8)          # N* = K0 at x0
  expect_equal(dd99_nd_fitness(xr, xr, n, p), 0, tolerance = 1e-12)

  d <- 1e-4
  g1 <- (dd99_nd_fitness(matrix(c(d, 0), 1), xr, n, p) -
           dd99_nd_fitness(matrix(c(-d, 0), 1), xr, n, p)) / (2 * d)
  g2 <- (dd99_nd_fitness(matrix(c(0, d), 1), xr, n, p) -
           dd99_nd_fitness(matrix(c(0, -d), 1), xr, n, p)) / (2 * d)
  expect_equal(c(g1, g2), c(0, 0), tolerance = 1e-6)          # gradient vanishes
  # branching (fitness minimum) in each direction since sigma_C < sigma_K
  curv1 <- (dd99_nd_fitness(matrix(c(d, 0), 1), xr, n, p) +
            dd99_nd_fitness(matrix(c(-d, 0), 1), xr, n, p)) / d^2
  expect_gt(curv1, 0)
})

test_that("DD99 2D multi-resident equilibrium zeroes resident fitness", {
  p <- list(r = 1, K0 = 500, x0 = c(0, 0), sigma_K = c(1, 1), sigma_C = c(0.5, 0.5))
  xr <- rbind(c(-1, 0), c(1, 0), c(0, 1))
  n <- dd99_nd_equilibrium(xr, p)
  expect_true(all(n > 0))
  expect_equal(dd99_nd_fitness(xr, xr, n, p), rep(0, 3), tolerance = 1e-8)
})

test_that("harness_dd99_nd validates parameter lengths", {
  expect_error(harness_dd99_nd(x0 = c(0, 0, 0), trait_names = c("x1", "x2")),
               "length")
})
