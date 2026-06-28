# Model harness (backend) abstraction for regnans
# ================================================
#
# A *harness* wires a `community` to a specific demographic model. The community
# pipeline (community_demography, community_selection_gradient,
# community_solve_singularity_1D, the assembler, ...) is model-agnostic: it only
# ever calls a small set of connector functions, each of which forwards to the
# harness's own implementation:
#
#   community_parameters()                    build model parameters
#   community_make_demography_runner()        closure: birth_rates -> offspring
#   community_demography_runner_cleanup()     write equilibrium state back
#   community_viable_bounds()                 viable trait region (empty community)
#   community_check_for_inviable_strategies() flag residents to drop
#   community_update_fitness_function()       build the invasion-fitness closure
#
# A harness is a list carrying a `$fns` table of these six functions plus any
# model state; its class is used for printing/dispatch. Two families ship:
#
#   harness_plant()    - the full `plant` SCM (the default; thin wrapper over the
#                        existing plant_community_* code). Versioned so different
#                        plant interfaces can be supported side by side.
#   harness_analytic() - fast toy models with analytic invasion fitness and
#                        equilibrium, used to develop and validate the assembly /
#                        attractor algorithms against known answers (issue #33).
#                        harness_bird() / harness_dd99() / harness_geritz98()
#                        are concrete instances.

# ---- connectors: the only entry points the pipeline uses --------------------

community_parameters <- function(community) {
  community$harness$fns$parameters(community)
}
community_make_demography_runner <- function(community) {
  community$harness$fns$make_demography_runner(community)
}
community_demography_runner_cleanup <- function(community, runner, converged = TRUE) {
  community$harness$fns$demography_runner_cleanup(community, runner, converged)
}
community_viable_bounds <- function(community) {
  community$harness$fns$viable_bounds(community)
}
community_check_for_inviable_strategies <- function(community) {
  community$harness$fns$check_for_inviable_strategies(community)
}
community_update_fitness_function <- function(community) {
  community$harness$fns$update_fitness_function(community)
}

# ---- plant harness ----------------------------------------------------------

##' The default model harness: the full `plant` SCM.
##'
##' Forwards each connector to the existing `plant_community_*` implementation,
##' so the working plant path is unchanged. The `version` argument is a hook for
##' supporting different plant interfaces side by side (currently informational).
##'
##' @title plant model harness
##' @param version plant-interface version label (informational for now)
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_plant <- function(version = "current") {
  h <- list(
    version = version,
    fns = list(
      parameters                    = plant_community_parameters,
      make_demography_runner        = plant_community_make_demography_runner,
      demography_runner_cleanup     = plant_community_demography_runner_cleanup,
      viable_bounds                 = plant_community_viable_bounds,
      check_for_inviable_strategies = plant_community_check_for_inviable_strategies,
      update_fitness_function       = plant_community_update_fitness_function
    )
  )
  class(h) <- c("harness_plant", "harness")
  h
}

# ---- generic analytic harness -----------------------------------------------

##' A harness for fast toy models with analytic invasion fitness.
##'
##' Implements all six connectors generically in terms of two model primitives:
##' a vectorised invasion-fitness function and a resident-equilibrium solver.
##' Concrete models (\code{harness_bird}, \code{harness_dd99},
##' \code{harness_geritz98}) supply these, typically backed by C++.
##'
##' Conventions:
##' \itemize{
##'   \item \code{community$birth_rate} stores resident abundance/density (the
##'     state the equilibrium solver returns), not a plant offspring rate.
##'   \item \code{fitness} returns LOG invasion fitness (resident value ~0).
##'   \item demography is solved by returning the analytic equilibrium directly;
##'     the standard \code{equilibrium_iteration} loop then converges in a couple
##'     of steps.
##' }
##'
##' @title Analytic (toy-model) harness
##' @param fitness function(x_mut, x_res, n_res, pars) -> numeric vector of log
##'   invasion fitness, one per mutant in \code{x_mut}.
##' @param equilibrium function(x_res, pars) -> numeric vector of resident
##'   equilibrium densities, one per resident.
##' @param pars named list of model parameters.
##' @param trait_names character vector naming the evolving trait(s).
##' @param label short model label (used in the harness class).
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_analytic <- function(fitness, equilibrium, pars, trait_names,
                             label = "analytic") {
  h <- list(
    pars        = pars,
    trait_names = trait_names,
    label       = label,
    # close pars into the model primitives so connectors call them with (x, ...)
    fitness     = function(x_mut, x_res, n_res) fitness(x_mut, x_res, n_res, pars),
    equilibrium = function(x_res) equilibrium(x_res, pars),
    fns = list(
      parameters                    = analytic_community_parameters,
      make_demography_runner        = analytic_community_make_demography_runner,
      demography_runner_cleanup     = analytic_community_demography_runner_cleanup,
      viable_bounds                 = analytic_community_viable_bounds,
      check_for_inviable_strategies = analytic_community_check_for_inviable_strategies,
      update_fitness_function       = analytic_community_update_fitness_function
    )
  )
  class(h) <- c(paste0("harness_", label), "harness_analytic", "harness")
  h
}

analytic_community_parameters <- function(community) {
  list(traits     = community$traits,
       birth_rate = community$birth_rate,
       pars       = community$harness$pars)
}

analytic_community_make_demography_runner <- function(community) {
  h <- community$harness
  last_offspring_production <- NULL
  history <- list()
  function(birth_rates) {
    x_res <- community$traits[, 1]
    out <- h$equilibrium(x_res)
    last_offspring_production <<- out
    history[[length(history) + 1L]] <<- list(`in` = birth_rates, out = out)
    out
  }
}

analytic_community_demography_runner_cleanup <- function(community, runner,
                                                         converged = TRUE) {
  e <- environment(runner)
  community$birth_rate <- e$last_offspring_production
  community$fitness_points <- NULL
  attr(community, "converged") <- converged
  attr(community, "progress") <- e$history
  community
}

analytic_community_update_fitness_function <- function(community) {
  h <- community$harness
  x_res <- community$traits[, 1]
  n_res <- community$birth_rate

  community$resident_fitness <-
    if (length(x_res) > 0L) h$fitness(x_res, x_res, n_res) else numeric()

  community$fitness_function <- function(x) {
    if (is.matrix(x)) x <- x[, 1]
    h$fitness(as.numeric(x), x_res, n_res)
  }
  community
}

analytic_community_viable_bounds <- function(community) {
  # Toy models use the bounds supplied to community_start(); we do not bracket a
  # fundamental-fitness region (these models are typically viable over the whole
  # supplied range). Bounds are already on the (empty) community.
  community
}

analytic_community_check_for_inviable_strategies <- function(community) {
  runner <- community_make_demography_runner(community)
  op <- runner(community$birth_rate)
  eps <- community$demography_control$equilibrium_extinct_birth_rate
  drop <- op < eps
  attr(op, "drop") <- drop
  op
}

# ---- concrete toy models ----------------------------------------------------

##' Migratory-bird arrival-time model (Johansson & Jonzen 2012; the simplified
##' analytic form of Brannstrom, Johansson & von Festenberg 2013, Games
##' 4:304-328, section 4).
##'
##' Trait = arrival time. Early arrival raises competitive ability
##' \code{C(x) = exp(-a x)} for a fixed number \code{K} of territories;
##' reproduction \code{R(x) = R0 exp(-(x - x_opt)^2 / (2 sigma^2))} peaks at the
##' seasonal optimum; survival \code{p}. The model has a single continuously
##' stable strategy (no branching) with the closed form
##' \deqn{x^* = x_{opt} - a\,\sigma^2,}
##' which arrives \emph{before} the population optimum (a tragedy of the commons).
##' This exact answer makes the model a benchmark for numerical
##' selection-gradient / singular-strategy solvers.
##'
##' @title Bird arrival-time harness
##' @param a competitive advantage of early arrival (>= 0)
##' @param x_opt seasonal optimum arrival time
##' @param sigma width of the benign season
##' @param R0 maximum reproductive output
##' @param K number of territories
##' @param p year-to-year survival, in (0, 1)
##' @param trait_name name of the evolving trait
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_bird <- function(a = 0.1, x_opt = 0, sigma = 1, R0 = 1, K = 1, p = 0.5,
                         trait_name = "x") {
  pars <- list(a = a, x_opt = x_opt, sigma = sigma, R0 = R0, K = K, p = p)
  harness_analytic(
    fitness     = bird_log_fitness,
    equilibrium = function(x_res, pars) bird_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "bird"
  )
}

##' Dieckmann & Doebeli 1999 competition model (Nature 400:354-357).
##'
##' Continuous-time logistic competition for a Gaussian resource: carrying
##' capacity \code{K(x) = K0 exp(-(x-x0)^2/(2 sigma_K^2))} and competition kernel
##' \code{C(d) = exp(-d^2/(2 sigma_C^2))}. Invasion fitness
##' \code{s(y) = r (1 - sum_i N_i C(y-x_i)/K(y))}. The singular strategy is
##' \code{x* = x0}; it is an evolutionary \emph{branching point} iff
##' \code{sigma_C < sigma_K} and an ESS iff \code{sigma_C > sigma_K} -- the
##' classic disruptive-selection oracle.
##'
##' @title Dieckmann & Doebeli 1999 harness
##' @param r intrinsic growth rate
##' @param K0 maximum carrying capacity
##' @param x0 trait optimising the resource (singular strategy)
##' @param sigma_K width of the resource/carrying-capacity kernel
##' @param sigma_C width of the competition kernel
##' @param trait_name name of the evolving trait
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_dd99 <- function(r = 1, K0 = 500, x0 = 0, sigma_K = 1, sigma_C = 0.4,
                         trait_name = "x") {
  pars <- list(r = r, K0 = K0, x0 = x0, sigma_K = sigma_K, sigma_C = sigma_C)
  harness_analytic(
    fitness     = dd99_fitness,
    equilibrium = function(x_res, pars) dd99_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "dd99"
  )
}

##' Geritz, Kisdi, Meszena & Metz 1998 soft-selection model (Evol. Ecol.
##' 12:35-57; the worked Levene example).
##'
##' \code{m} patches with optima \code{mu} and capacities \code{K}; Gaussian
##' within-patch survival of width \code{sigma}. Invasion fitness reduces, for a
##' single resident, to \code{S(y) = log sum_j c_j f_j(y)/f_j(x)} with
##' \code{c_j = K_j/sum K}. For the symmetric three-patch default
##' (\code{mu = (-d, 0, d)}, equal \code{K}) the singular strategy is
##' \code{x* = 0}: convergence stable for all \code{d/sigma}, and an evolutionary
##' branching point iff \code{d/sigma > sqrt(3/2) ~= 1.2247} (a CSS below that).
##'
##' @title Geritz et al. 1998 soft-selection harness
##' @param d patch-optimum spacing for the symmetric default (optima -d, 0, d)
##' @param sigma within-patch survival width
##' @param mu patch optima (overrides the symmetric default built from \code{d})
##' @param K patch capacities (defaults to equal capacities)
##' @param trait_name name of the evolving trait
##' @return a `harness` object
##' @author Daniel Falster
##' @export
harness_geritz98 <- function(d = 1.5, sigma = 1, mu = c(-d, 0, d),
                             K = rep(1, length(mu)), trait_name = "x") {
  pars <- list(sigma = sigma, mu = mu, K = K)
  harness_analytic(
    fitness     = geritz_log_fitness,
    equilibrium = function(x_res, pars) geritz_equilibrium(x_res, pars),
    pars        = pars,
    trait_names = trait_name,
    label       = "geritz98"
  )
}

##' @export
print.harness <- function(x, ...) {
  cat(sprintf("<harness: %s>\n", paste(class(x)[-length(class(x))], collapse = ", ")))
  if (!is.null(x$pars)) {
    cat("  parameters:", paste(names(x$pars), unlist(x$pars), sep = "=",
                               collapse = ", "), "\n")
  }
  invisible(x)
}
