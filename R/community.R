##' Initialise a community object
##'
##' Used to store full description of community. Is a list
##' with elements parameters, bounds, birth_rate_initial,
##' trait_names, traits, birth_rate, fitness_control.
##'
##' @title Initialise a community object
##' @param parameters A \code{parameters} object, as specified
##' in \code{plant}..
##' @param bounds A set of bounds, as specified in \code{plant}.
##' @param birth_rate_initial A vector of birth rates.
##' @param hyperpar A plant hyperparameter function to be used when calling \code{strategy_list}
##' @param fitness_control List of parameters controlling
##' how approximate fitness landscapes are generated. See function
##' \code{fitness_control} for an example.
##' @return An \code{community} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
## TODO: Put birth_rate_initial into parameters and set up
## appropriately?  We use 1 in a couple of places, no?
community_start <- function(bounds,
                            birth_rate_initial = 1e-3,
                            model_support = NULL,
                            fitness_control = NULL) {

  if (is.character(bounds)) {
    bounds <-  bounds_infinite(bounds)
  }
  ## TODO: Check parameters is empty.
  ret <- list(bounds = check_bounds(bounds),
              birth_rate_initial = birth_rate_initial,
              model_support = model_support
              )
  ret$trait_names <- rownames(bounds)
  ret$traits <- trait_matrix(numeric(0), ret$trait_names)
  ret$birth_rate <- numeric(0)
  ret$fitness_control <- fitness_control
  class(ret) <- "community"
  ret
}

community_add <- function(obj, traits, birth_rate=NULL) {

  if (is.null(birth_rate)) {
    birth_rate <- obj$birth_rate_initial
  }
  if (!is.matrix(traits)) {
    stop("traits must be a matrix")
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

    obj <- community_reset(obj)
  }
  obj
}

community_drop <- function(community, which) {
  n_spp <- length(community)
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
    community$traits <- community$traits[keep,,drop=FALSE]
    community$birth_rate <- community$birth_rate[keep]
    community <- community_reset(community)
  }

  if (is.null(community$fitness_function)) {
    community <- community %>% community_run()
  }

  community
}

community_reset <- function(community) {
  community$resident_fitness <- NULL
  community$fitness_points <- NULL
  community$fitness_slopes <- NULL
  community$fitness_function <- NULL
  
  community$model_support$node_schedule_times <- NULL

  community
}


##' Helper function to create trait matrices suitable for
##' \code{\link{strategy_list}}.
##'
##' @title Create trait matrix
##' @param x Values
##' @param trait_name Name of a single trait
##' @export
##' @author Rich FitzJohn
trait_matrix <- function(x, trait_name) {
  m <- matrix(x, ncol = length(trait_name))
  colnames(m) <- trait_name
  m
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
