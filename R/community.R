## TODO: Put seed_rain_initial into parameters and set up
## appropriately?  We use 1 in a couple of places, no?
##' @export
community <- function(parameters, bounds, seed_rain_initial=1e-3,
                      fitness_approximate_control=NULL) {
  if (is.character(bounds)) {
    bounds <- bounds_infinite(bounds)
  }
  ## TODO: Check parameters is empty.
  ret <- list(parameters=validate(parameters),
              bounds=check_bounds(bounds),
              seed_rain_initial=seed_rain_initial)
  ret$trait_names <- rownames(bounds)
  ret$traits <- trait_matrix(numeric(0), ret$trait_names)
  ret$seed_rain <- numeric(0)
  ret$fitness_approximate_control <- fitness_approximate_control
  class(ret) <- "community"
  ret
}

make_community <- function(traits, seed_rain, parameters,
                           bounds=colnames(traits), ...) {
  community_add(community(parameters, bounds, ...),
                traits, seed_rain)
}

community_parameters <- function(obj) {
  p <- obj$parameters
  p$strategies <- strategy_list(obj$traits, p)
  p$seed_rain <- obj$seed_rain

  if (!is.null(obj$cohort_schedule_times)) {
    p$cohort_schedule_times <- obj$cohort_schedule_times
  }
  if (!is.null(obj$cohort_schedule_ode_times)) {
    p$cohort_schedule_ode_times <- obj$cohort_schedule_ode_times
  }

  p
}

community_viable_bounds <- function(obj) {
  if (length(obj) > 0) {
    stop("You don't want to run this on an existing community")
  }
  obj$bounds <- viable_fitness(obj$bounds, community_parameters(obj))
  obj
}

community_add <- function(obj, traits, seed_rain=NULL) {
  if (is.null(seed_rain)) {
    seed_rain <- obj$seed_rain_initial
  }
  if (!is.matrix(traits)) {
    stop("traits must be a matrix") # sensible?
  }
  if (length(seed_rain) == 1) {
    seed_rain <- rep_len(seed_rain, nrow(traits))
  }
  if (length(seed_rain) != nrow(traits)) {
    stop("Incompatible length seed rain")
  }
  if (ncol(traits) != ncol(obj$traits)) {
    stop("Incorrect size trait matrix")
  }
  if (!is.null(colnames(traits)) &&
      !identical(colnames(traits), obj$trait_names)) {
    stop("Incorrect traits")
  }
  if (nrow(traits) > 0L) {
    obj$traits <- rbind(obj$traits, traits)
    obj$seed_rain <- c(obj$seed_rain, seed_rain)
    ## Need to deal with cohort times here.  I think what we should do
    ## here is to retain the times as best we can?  For now though,
    ## we'll just nuke the times.
    obj <- community_clear_times(obj)
  }
  obj
}

community_drop <- function(obj, which) {
  n_spp <- length(obj)
  if (is.logical(which)) {
    if (length(which) != n_spp) {
      stop(sprintf("Invalid length: expected %d, recieved %d",
                   n_spp, length(which)))
    }
    keep <- !which
  } else if (is.numeric(which) || is.integer(which)) {
    which <- as.integer(which)
    if (any(which < 1 || which > n_spp)) {
      stop("Invalid indicies")
    }
    keep <- !(which %in% seq_len(n_spp))
  } else {
    stop("Invalid index")
  }
  if (!all(keep)) {
    obj$traits <- obj$traits[keep,,drop=FALSE]
    obj$seed_rain <- obj$seed_rain[keep]
    obj <- community_clear_times(obj)
  }
  obj
}

community_clear_times <- function(obj) {
  obj$cohort_schedule_ode_times <- NULL
  obj$cohort_schedule_times <- NULL
  obj$fitness_approximate_points <- NULL
  obj$fitness_approximate_slopes <- NULL
  obj
}

community_run <- function(obj) {
  if (length(obj) > 0L) {
    p <- build_schedule(community_parameters(obj))
    obj$seed_rain <- attr(p, "seed_rain_out")
    obj$cohort_schedule_times <- p$cohort_schedule_times
    obj$cohort_schedule_ode_times <- p$cohort_schedule_ode_times
    obj$fitness_approximate_points <- NULL
  }
  obj
}

community_run_to_equilibrium <- function(obj) {
  if (length(obj) > 0L) {
    p <- equilibrium_seed_rain(community_parameters(obj))
    obj$seed_rain <- attr(p, "seed_rain_out")
    obj$cohort_schedule_times <- p$cohort_schedule_times
    obj$cohort_schedule_ode_times <- p$cohort_schedule_ode_times
    obj$fitness_approximate_points <- NULL
  }
  obj
}

community_fitness <- function(obj, traits) {
  community_make_fitness(obj)(traits)
}

##' @export
length.community <- function(x) {
  nrow(x$traits)
}
