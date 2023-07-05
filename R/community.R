##' Initialise a community object
##'
##' Used to store full description of community. Is a list
##' with elements parameters, bounds, birth_rate_initial,
##' trait_names, traits, birth_rate, fitness_approximate_control.
##'
##' @title Initialise a community object
##' @param parameters A \code{parameters} object, as specified
##' in \code{plant}..
##' @param bounds A set of bounds, as specified in \code{plant}.
##' @param birth_rate_initial A vector of birth rates.
##' @param hyperpar A plant hyperparameter function to be used when calling \code{strategy_list}
##' @param fitness_approximate_control List of parameters controlling
##' how approximate fitness landscapes are generated. See function
##' \code{fitness_approximate_control} for an example.
##' @return An \code{community} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
## TODO: Put birth_rate_initial into parameters and set up
## appropriately?  We use 1 in a couple of places, no?
community_start <- function(parameters, bounds,
                            birth_rate_initial = 1e-3,
                            hyperpar = NULL,
                            fitness_approximate_control = NULL) {

  if (is.character(bounds)) {
    bounds <- bounds_infinite(bounds)
  }
  ## TODO: Check parameters is empty.
  ret <- list(parameters=validate(parameters),
              bounds=check_bounds(bounds),
              birth_rate_initial=birth_rate_initial, 
              hyperpar = hyperpar
              )
  ret$trait_names <- rownames(bounds)
  ret$traits <- trait_matrix(numeric(0), ret$trait_names)
  ret$birth_rate <- numeric(0)
  ret$fitness_approximate_control <- fitness_approximate_control
  class(ret) <- "community"
  ret
}

make_community <- function(traits, birth_rate, parameters,
                           bounds=colnames(traits), ...) {
  community_add(community(parameters, bounds, ...),
                traits, birth_rate)
}

community_parameters <- function(obj) {
  
  p <- obj$parameters

  if(is.null(obj$hyperpar))
    hyperpar <- param_hyperpar(p)
  else 
    hyperpar <- obj$hyperpar
  
  p$strategies <- strategy_list(obj$traits, p, birth_rate_list = obj$birth_rate, hyperpar = hyperpar)

  if (!is.null(obj$node_schedule_times)) {
    p$node_schedule_times <- obj$node_schedule_times
  }
  if (!is.null(obj$node_schedule_ode_times)) {
    p$node_schedule_ode_times <- obj$node_schedule_ode_times
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

community_add <- function(obj, traits, birth_rate=NULL) {
  if (is.null(birth_rate)) {
    birth_rate <- obj$birth_rate_initial
  }
  if (!is.matrix(traits)) {
    stop("traits must be a matrix") # sensible?
  }
  if (length(birth_rate) == 1) {
    birth_rate <- rep_len(birth_rate, nrow(traits))
  }
  if (length(birth_rate) != nrow(traits)) {
    stop("Incompatible length birth rate")
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
    obj$birth_rate <- c(obj$birth_rate, birth_rate)
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
    obj$birth_rate <- obj$birth_rate[keep]
    obj <- community_clear_times(obj)
  }
  obj
}

community_clear_times <- function(obj) {
  obj$node_schedule_ode_times <- NULL
  obj$node_schedule_times <- NULL
  obj$fitness_approximate_points <- NULL
  obj$fitness_approximate_slopes <- NULL
  obj
}

community_run <- function(obj) {
  if (length(obj) > 0L) {
    p <- build_schedule(community_parameters(obj), ctrl = plant_control())
    obj$birth_rate <- attr(p, "offspring_production")
    obj$node_schedule_times <- p$node_schedule_times
    obj$node_schedule_ode_times <- p$node_schedule_ode_times
    obj$fitness_approximate_points <- NULL
  }
  obj
}

community_run_to_equilibrium <- function(obj) {
  if (length(obj) > 0L) {
    p <- equilibrium_birth_rate(community_parameters(obj), 
            ctrl = plant_control())

    obj$birth_rate <- attr(p, "offspring_production")
    obj$node_schedule_times <- p$node_schedule_times
    obj$node_schedule_ode_times <- p$node_schedule_ode_times
    obj$fitness_approximate_points <- NULL
  }
  obj
}

community_fitness <- function(obj, traits) {
  community_make_fitness(obj)(traits)
}

##' Returns number of types in community
##'
##' @title Returns number of types in \code{community}
##' @param x A \code{community} object.
##' @return number of rows in trait matrix.
##' @author Rich FitzJohn, Daniel Falster
##' @export
length.community <- function(x) {
  nrow(x$traits)
}
