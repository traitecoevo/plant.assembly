##' Fitness of trait value(s) evaluated against a community.
##'
##' Returns the invasion fitness (log net reproduction ratio) of one or more
##' mutant trait values in the community's environment. For an *empty* community
##' this is the fundamental fitness (growth with no competitors).
##'
##' Reimplemented on the community machinery: it evaluates
##' \code{community$fitness_function} (built by
##' \code{plant_community_update_fitness_function}), replacing plant's removed
##' \code{fundamental_fitness()} / \code{fitness_landscape()}.
##'
##' @title Fitness of trait values against a community
##' @param community A \code{community} object.
##' @param values A vector (1D) or matrix (multi-trait) of trait values, as for
##'   \code{\link{trait_matrix}}.
##' @return A numeric vector of fitnesses, one per trait value.
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_growth_rate <- function(community, values) {
  if (is.null(community$fitness_function)) {
    community <- community_update_fitness_function(community)
  }
  community$fitness_function(values)
}

##' Find the trait value of maximum fitness within some bounds.
##'
##' Searches \code{bounds} for the trait value(s) that maximise invasion fitness
##' into the community's environment (the fundamental niche peak for an empty
##' community). Uses \code{stats::optimise} in 1D and \code{stats::optim}
##' (L-BFGS-B) in higher dimensions, operating on \code{community$fitness_function}.
##'
##' @title Find point of maximum fitness within some range
##' @param community A \code{community} object.
##' @param bounds Bounds matrix (\code{lower}/\code{upper} per trait). Defaults
##'   to \code{community$bounds}.
##' @param log_scale Is the trait naturally on a log scale? If so the search is
##'   done in log space, which is usually much better behaved.
##' @param tol Tolerance passed to the optimiser.
##' @return The maximising trait value(s), named by trait, with the achieved
##'   fitness in attribute \code{"fitness"}.
##' @importFrom stats optimise optim
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_fitness <- function(community, bounds = NULL, log_scale = TRUE, tol = 1e-3) {
  if (is.null(bounds)) {
    bounds <- community$bounds
  }
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)

  if (is.null(community$fitness_function)) {
    community <- community_update_fitness_function(community)
  }
  fitness <- community$fitness_function

  if (log_scale) {
    bounds[bounds[, 1] == -Inf, 1] <- 0
    bounds <- log(bounds)
    f <- function(x) fitness(exp(x))
  } else {
    f <- function(x) fitness(x)
  }

  ret <- solve_max_worker(bounds, f, tol = tol)
  loc <- as.numeric(ret)
  if (log_scale) {
    loc <- exp(loc)
  }
  structure(loc, names = traits, fitness = attr(ret, "fitness"))
}

## Numeric maximiser used by max_fitness(): `f` takes trait values on the
## (already transformed) search scale and returns a scalar fitness; `bounds` are
## the search bounds on that same scale. Returns the maximising location with
## the achieved value in attribute "fitness".
solve_max_worker <- function(bounds, f, tol = 1e-3) {
  if (nrow(bounds) == 1L) {
    if (!all(is.finite(bounds))) {
      stop("Starting value did not have finite fitness; finite bounds required")
    }
    ## suppressWarnings: optimise warns "NA/Inf replaced by maximum positive
    ## value" for inviable trait values, which is the behaviour we want.
    out <- suppressWarnings(optimise(f, interval = bounds, maximum = TRUE, tol = tol))
    structure(out$maximum, fitness = out$objective)
  } else {
    ## Not very well tested, and the tolerance is not useful:
    out <- optim(rowMeans(bounds), f, method = "L-BFGS-B",
                 lower = bounds[, "lower"], upper = bounds[, "upper"],
                 control = list(fnscale = -1, factr = 1e10))
    structure(out$par, fitness = out$value)
  }
}
