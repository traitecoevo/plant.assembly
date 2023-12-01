assembler_births <- function(obj) {
  to_add <- assembler_new_types(obj)

  ## Now we add things:
  if (nrow(to_add) > 0) {
    plant_log_births_new_types(to_add)
    i <- should_move(obj, to_add)
    if (length(i) == 1L) {
      obj <- assembler_births_try_move(obj, to_add, i)
    } else {
      rain <- initial_birth_rate(attr(to_add, "fitness"), obj)
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
    obj$community$done <- TRUE
  }
  obj
}

assembler_births_try_move <- function(obj, to_add, i) {
  community <- obj$community
  plant_log_births(paste0("Trying to move resident ", i))
  x_prev <- community$traits[i, , drop=FALSE]

  test <- community
  test$traits[i, ] <- to_add
  test <- community_run_to_equilibrium(test)

  ## Now, try adding the previous case back in.  First, check
  ## the fitness:
  w_prev <- test$fitness_function(x_prev)
  if (w_prev > 0) {
    ## Here, there is no need to run anything: the original
    ## resident has positive fitness and will increase
    plant_log_births("Adding original resident back")
    community <- community_add(test, x_prev, initial_birth_rate(w_prev, obj))
  } else {
    plant_log_births("Move was successful")
    community <- test
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

## First order birth rate approximation based on fitness.
initial_birth_rate <- function(fitness, obj) {
  if (is.null(fitness) || obj$control$run_type == "single") {
    birth_rate <- NULL
  } else {
    birth_rate <- pmin(exp(fitness), obj$control$max_birth_rate_initial)
    birth_rate <- pmax(birth_rate, obj$control$min_birth_rate_initial)
  }
  birth_rate
}

plant_log_births_new_types <- function(to_add) {
  str <- format_community_state(to_add, "-", sprintf("\t%s\t", c("i", seq_along(to_add))))
  plant_log_births(paste(c("Proposed new type(s):", str), collapse="\n"),
                   to_add=to_add)
}
