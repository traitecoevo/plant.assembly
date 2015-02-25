assembler_births <- function(obj) {
  to_add <- assembler_new_types(obj)

  ## Now we add things:
  if (nrow(to_add) > 0) {
    message_new_types(to_add)
    i <- should_move(obj, to_add)
    if (length(i) == 1L) {
      obj <- assembler_births_try_move(obj, to_add, i)
    } else {
      rain <- initial_seed_rain(attr(to_add, "fitness"), obj)
      obj$community <- community_add(obj$community, to_add, rain)
    }
  }

  ## This is a bit ugly: we only really want to respond to no births
  ## with being done if we are running to equilibrium.  If we're using
  ## the simple step, there is no way of finishing here unless we look
  ## to a limited amount of change in the community since last step.
  if (has_attr(to_add, "done")) {
    obj$done <- (obj$control$run_type == "to_equilibrium" &&
                 attr(to_add, "done"))
  }
  obj
}

assembler_births_try_move <- function(obj, to_add, i) {
  community <- obj$community
  message("assembler[births]> Trying to move resident ", i)
  test <- community
  prev <- test$traits
  test$traits[i,] <- to_add
  test <- community_run_to_equilibrium(test)

  ## Now, try adding the previous case back in.  First, check
  ## the fitness:
  fitness <- community_make_fitness(test)
  ## TODO: extract times from this.
  w <- fitness(test$traits[i, , drop=FALSE])
  if (w > 0) {
    ## Here, there is no need to run anything: the original
    ## resident has positive fitness and will increase
    message("assembler[births]> Adding original resident back")
    community <-
      community_add(test, test$traits[i, , drop=FALSE],
                    initial_seed_rain(w, obj))
  } else {
    message("assembler[births]> Move was successful")
    community <- test
    ## Most times are fine, but the approximate points are going to be
    ## broken.
    community$fitness_approximate_points <- NULL
  }
  obj$community <- community
  obj
}

assembler_new_types <- function(obj) {
  community_new_types(obj$community, obj$control)
}

community_new_types <- function(community, control) {
  if (control$birth_type == "maximum") {
    community_new_types_maximum_fitness(community, control)
  } else if (control$birth_type == "stochastic") {
    community_new_types_stochastic(community, control)
  } else {
    stop("Unknown birth type ", control$birth_type)
  }
}

## First order seed rain approximation based on fitness.
initial_seed_rain <- function(fitness, obj) {
  if (is.null(fitness) || obj$control$run_type == "single") {
    seed_rain <- NULL
  } else {
    ## This magic number is not great, but needed to prevent
    ## suggesting adding 1e8 as the seed rain (can happen!).  This
    ## should be "quite large", but obviously that number depends on
    ## the situation.
    max_seed_rain_initial <- 500
    min_seed_rain_initial <-
      obj$community$parameters$control$equilibrium_extinct_seed_rain
    seed_rain <- pmin(exp(fitness), max_seed_rain_initial)
    seed_rain <- pmax(seed_rain, min_seed_rain_initial)
  }
  seed_rain
}
