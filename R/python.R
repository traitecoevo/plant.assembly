##' This is a small helper function for calling tree from Python.  As
##' we develop the ideas better, this will change.  The current verison
##' allows computing of the fitness landscape for a set of parameters,
##' similar to the test case we've been talking about.
##' @title Wrapper code for Python
##' @param time_disturbance The mean time to disturbance
##' @param slope Indicates strength of a leaf economic spectrum tradeoff
##' @param lma Vector of species traits
##' @param seed_rain Vector of species abundances
##' @param equilibrium Logical, inidicating if we should compute
##' equilibrium seed rain (the default) or not.
##' @param equilibrium_nsteps Maximum number of steps to run the
##' equilibrium rain for
##' @param verbose Be noisy during calculation of equilibrium seed
##' rain.
##' @export
##' @rdname for_python
for_python_equilibrium <- function(time_disturbance, slope,
                                   lma, seed_rain,
                                   equilibrium_nsteps=20L, verbose=FALSE) {
  if (length(seed_rain) != length(lma)) {
    stop("seed_rain must be same length as lma")
  }

  sys0 <- python_base_community(time_disturbance, slope, verbose,
                                equilibrium_nsteps)
  sys0 <- community_add(sys0, trait_matrix(lma, "lma"), seed_rain)
  sys0 <- community_run_to_equilibrium(sys0)
  sys0 <- community_prepare_approximate_fitness(sys0)
  sys0
}

## Based on successional_diversity:analysis/R/assembly.R:run_model()
##' @param nsteps Number of evolutionary steps to run for
##' @export
##' @rdname for_python
for_python_evolve <- function(time_disturbance, slope, nsteps,
                              filename=NULL, lma=NULL, seed_rain=NULL,
                              verbose=FALSE) {
  equilibrium_nsteps <- 20L
  sys0 <- python_base_community(time_disturbance, slope,
                                verbose, equilibrium_nsteps)
  max_bounds <- sys0$bounds
  obj <- assembler(sys0, list(birth_type="maximum", birth_move_tol=1),
                   filename=filename)
  if (!is.null(lma)) {
    obj <- assembler_set_traits(obj, trait_matrix(lma, "lma"), seed_rain)
  }
  obj <- assembler_run(obj, nsteps)
  if (is.null(obj$community$bounds)) {
    obj$community$bounds <- max_bounds
  }

  obj$community
}

##' @export
##' @param filename Optional filename to save results in (as for
##'   \code{\link{assembler}}).
##' @rdname for_python
for_python_simulate <- function(time_disturbance, slope, traits,
                                nsteps=100, filename=NULL) {
  sys <- for_python_evolve(time_disturbance, slope, nsteps, filename)
  for_python_fitness(sys, traits)
}

python_base_community <- function(time_disturbance, slope,
                                  verbose=FALSE,
                                  equilibrium_nsteps=20) {
  p0 <- ebt_base_parameters()
  if (verbose) {
    p0$control <- equilibrium_verbose(p0$control)
  } else {
    p0$control <- equilibrium_quiet(p0$control)
  }
  p0$control$equilibrium_nsteps  <- equilibrium_nsteps
  p0$control$equilibrium_solver_name <- "hybrid"

  p0$strategy_default$c_r1 <- 0.5
  p0$strategy_default$c_r2 <- 0
  p0$hyperpar <- make_FFW16_hyperpar(B_kl2=slope)

  p0$disturbance_mean_interval <- time_disturbance

  control <- list(type="gp", cost=grail::cost_mean_2sd)

  sys0 <- community(p0, bounds(lma=10^c(-4, 2)),
                    fitness_approximate_control=control)

  sys0
}

##' @export
##' @param sys Result of running \code{for_python} or
##' \code{for_python_evolve}
##' @rdname for_python
for_python_fitness <- function(sys, lma) {
  community_make_fitness(sys)(trait_matrix(lma, "lma"))
}
##' @export
##' @rdname for_python
for_python_fitness_approximate <- function(sys, lma) {
  community_fitness_approximate(sys)(trait_matrix(lma, "lma"))
}
##' @export
##' @rdname for_python
for_python_fitness_prepare <- function(sys) {
  community_prepare_approximate_fitness(sys)
}
