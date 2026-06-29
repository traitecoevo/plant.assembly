# Plant SCM smoke tests (issue #27).
#
# The model-agnostic pipeline (demography, selection gradient, fitness
# landscape, attractor solving, assembly, the alternative equilibrium solvers,
# util_nlsolve, inviable-strategy checks) is exercised fast on the DD99 toy
# harness elsewhere. This file keeps a *minimal* set of genuine plant-SCM runs
# as the integration anchor: it is the only place that confirms the real
# physiological model still drives the regnans pipeline, and it pins the
# reference values for the installed plant (2.0.0.9001).
#
# These are the slow tests in the suite (a few seconds each, max_patch_lifetime
# = 30 via helper-assembly.R). Keep the count small.

test_that("empty community: fundamental niche and fitness peak (SCM)", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support())

  # fundamental niche: the positive-fitness interval, strictly inside the bounds
  vb <- community_viable_fitness_1D(comm)
  expect_equal(rownames(vb), "lma")
  expect_equal(colnames(vb), c("lower", "upper"))
  expect_lt(vb[, "lower"], vb[, "upper"])
  # reference values for the installed plant (2.0.0.9001)
  expect_equal(unname(vb[, "lower"]), 0.04045, tolerance = 1e-3)
  expect_equal(unname(vb[, "upper"]), 0.79930, tolerance = 1e-3)

  ff <- plant_community_update_fitness_function(comm)$fitness_function
  expect_gt(ff(exp(mean(log(vb)))), 1)          # clearly positive inside
  expect_lt(abs(ff(vb[, "lower"])), 0.01)       # ~0 at the edges
  expect_lt(abs(ff(vb[, "upper"])), 0.01)

  # max_fitness / max_growth_rate against the same fitness function
  mx <- max_fitness(comm)
  expect_equal(names(mx), "lma")
  expect_gt(as.numeric(mx), 0.01)
  expect_lt(as.numeric(mx), 2)
  g <- max_growth_rate(comm, c(0.05, 0.0825, 0.2))
  expect_equal(g, ff(c(0.05, 0.0825, 0.2)), tolerance = 1e-8)
  grid <- max_growth_rate(comm, seq_log_range(c(0.01, 2), 7))
  expect_gte(attr(mx, "fitness"), max(grid) - 1e-3)
})

test_that("single resident solves to equilibrium, gradient and inviable check (SCM)", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography()

  expect_true(attr(comm, "converged"))
  # equilibrium birth rate, reference for plant 2.0.0.9001 (the genuinely
  # iterative equilibrium_iteration fixed point, robust to the starting rate)
  expect_equal(as.numeric(comm$birth_rate), 0.068455, tolerance = 1e-3)
  expect_equal(comm$resident_fitness, 0, tolerance = 1e-4)

  # selection gradient is finite at the (non-singular) resident
  comm <- community_selection_gradient(comm)
  expect_length(comm$selection_gradient, 1L)
  expect_true(is.finite(comm$selection_gradient))

  # plant_community_check_for_inviable_strategies: the viable resident is kept
  op <- community_check_for_inviable_strategies(comm)
  expect_false(any(attr(op, "drop")))
})

test_that("community_solve_singularity_1D finds the 1D attractor (SCM)", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_solve_singularity_1D(bounds = c(0.13, 0.16), tol = 1e-3)

  expect_length(comm, 1L)
  x <- as.numeric(comm$traits[, "lma"])
  expect_equal(x, 0.1417, tolerance = 1e-2)     # reference for plant 2.0.0.9001
  # the gradient is small relative to its scale (root tol is on the trait)
  expect_true(is.finite(comm$selection_gradient))
  expect_lt(abs(comm$selection_gradient), 50)
})
