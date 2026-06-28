# Tests for the community interface (community_start/add/drop, demography,
# selection gradient). These target regnans's own API, unlike the older
# test files which call plant's (now removed) internal functions directly.

# `assembly_model_support()` is defined in helper-assembly.R.

# ---- fast unit tests (no SCM runs) ----------------------------------------

test_that("trait_matrix builds a named single-trait matrix", {
  m <- trait_matrix(c(0.1, 0.2), "lma")
  expect_true(is.matrix(m))
  expect_equal(dim(m), c(2L, 1L))
  expect_equal(colnames(m), "lma")
  expect_equal(as.numeric(m), c(0.1, 0.2))
})

test_that("bounds() produces lower/upper rows named by trait", {
  b <- bounds(lma = c(0.01, 2))
  expect_equal(rownames(b), "lma")
  expect_equal(unname(b["lma", ]), c(0.01, 2))
})

test_that("demographic_step_control has expected defaults and rejects unknowns", {
  ctrl <- demographic_step_control()
  expect_equal(ctrl$equilibrium_solver_name, "equilibrium_iteration")
  for (nm in c("equilibrium_eps", "equilibrium_nsteps",
               "equilibrium_large_birth_rate_change",
               "equilibrium_extinct_birth_rate",
               "equilibrium_min_offspring_arriving")) {
    expect_true(nm %in% names(ctrl), info = nm)
  }
  expect_equal(demographic_step_control(list(equilibrium_eps = 1e-3))$equilibrium_eps,
               1e-3)
  expect_error(demographic_step_control(list(not_a_param = 1)),
               "Unknown control parameters")
})

test_that("community_start creates an empty community", {
  comm <- community_start(bounds(lma = c(0.01, 2)))
  expect_s3_class(comm, "community")
  expect_equal(comm$trait_names, "lma")
  expect_equal(length(comm), 0L)
  expect_equal(nrow(comm$traits), 0L)
})

test_that("community_add appends species and validates input", {
  comm <- community_start(bounds(lma = c(0.01, 2)))
  comm <- community_add(comm, trait_matrix(c(0.05, 0.1), "lma"), birth_rate = 5)
  expect_equal(length(comm), 2L)
  expect_equal(comm$birth_rate, c(5, 5))                 # recycled scalar
  expect_error(community_add(comm, c(0.05, 0.1)), "must be a matrix")
  expect_error(community_add(comm, trait_matrix(0.05, "lma"), birth_rate = c(1, 2)),
               "Incompatible length")
  expect_error(community_add(comm, trait_matrix(c(0.05, 0.1), c("lma", "x"))),
               "Incorrect size")
})

test_that("community_drop validates indices (no solve needed)", {
  # Index validation happens before any demography re-solve, so no model_support.
  comm <- community_start(bounds(lma = c(0.01, 2)))
  comm <- community_add(comm, trait_matrix(c(0.05, 0.1, 0.2), "lma"),
                        birth_rate = c(1, 2, 3))
  expect_error(community_drop(comm, 99), "Invalid indicies")
  expect_error(community_drop(comm, c(TRUE, FALSE)), "Invalid length")
})

test_that("community_drop keeps the correct species (regression: numeric index)", {
  # community_drop resets and re-solves the survivors, so model_support is
  # required. Traits are stable under the re-solve (only birth rates change), so
  # we assert on which traits survive.
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_add(trait_matrix(c(0.05, 0.2), "lma"), birth_rate = c(50, 50))

  d <- community_drop(comm, 2)            # drop the 2nd species
  expect_equal(length(d), 1L)
  expect_equal(as.numeric(d$traits[, "lma"]), 0.05)

  d2 <- community_drop(comm, c(FALSE, TRUE))  # logical mask, drop the 2nd
  expect_equal(as.numeric(d2$traits[, "lma"]), 0.05)
})

# ---- regression: schedule defaults track max_patch_lifetime ----------------

test_that("plant_default_assembly_pars regenerates the node schedule default", {
  # Lowering max_patch_lifetime must shrink node_schedule_times_default to match,
  # otherwise the SCM errors with "time_max must be greater than current time"
  # when the equilibrium runner resets the schedule.
  p <- plant_default_assembly_pars(max_patch_lifetime = 30)
  expect_equal(p$max_patch_lifetime, 30)
  expect_lte(max(p$node_schedule_times_default), 30)
})

# ---- integration tests (run the SCM) ---------------------------------------

test_that("community_demography on an empty community gives no resident fitness", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support())
  comm <- community_demography(comm)
  expect_length(comm$resident_fitness, 0L)
  expect_true(is.function(comm$fitness_function))
})

test_that("community_demography solves a single resident to equilibrium", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography()

  expect_true(attr(comm, "converged"))
  expect_length(comm$birth_rate, 1L)
  # equilibrium birth rate (reference from current plant 2.0.0.9001; robust to
  # the starting birth rate — converges to the same fixed point from 20/200/1000)
  expect_equal(comm$birth_rate, 0.068455, tolerance = 1e-3)
  # at equilibrium the resident's invasion fitness is ~0
  expect_equal(comm$resident_fitness, 0, tolerance = 1e-4)
})

test_that("community_selection_gradient returns a finite gradient at the resident", {
  comm <- community_start(bounds(lma = c(0.01, 2)),
                          model_support = assembly_model_support()) |>
    community_add(trait_matrix(0.0825, "lma"), birth_rate = 200) |>
    community_demography() |>
    community_selection_gradient()

  expect_length(comm$selection_gradient, 1L)
  expect_true(is.finite(comm$selection_gradient))
})
