assembler <- function(community, control=NULL, filename=NULL, prev=NULL) {
  control <- assembler_control(control)

  ret <- list(births=make_births(control),
              deaths=make_deaths(control),
              history=list(),
              done=FALSE,
              control=control)

  if (is.null(prev)) {
    ret <- assembler_initialise(ret, community)
  } else {
    ret <- assembler_restore(ret, prev)
  }

  class(ret) <- "assembler"
  ret
}

assembler_control <- function(control=NULL) {
  defaults <- list(run_type="to_equilibrium",
                   birth_type="maximum",
                   birth_move_tol=0,
                   compute_viable_fitness=TRUE,
                   ## births_stochastic_naive:
                   n_mutants=1L,
                   n_immigrants=1L,
                   check_positive=TRUE,
                   check_inviable=TRUE,
                   ## + vcv, which has no default.
                   ## deaths_stochastic_naive:
                   dead_seed_rain=1e-3,
                   ## births_maximum_fitness:
                   eps_too_close=1e-3)
  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  stop("Unknown control parameters ", paste(extra, collapse=", "))
  ret <- modifyList(defaults, control)
  if (ret$birth_type == "maximum" && ret$run_type != "to_equilibrium") {
    stop("Must use 'to_equilibrium' run type with maximum fitness births")
  }
  ret
}

make_births <- function(control) {
  switch(control$birth_type,
         stochastic=make_births_stochastic_naive(control),
         maximum=make_births_maximum_fitness(control),
         stop("Unknown births type ", control$type))
}

make_deaths <- function(control) {
  make_deaths_stochastic_naive(control)
}

assembler_initialise <- function(obj, community) {
  if (!inherits(community, "community")) {
    stop("Expecting a community object")
  }
  message("assembler> Starting empty assembler")
  obj$community <- community
  if (isTRUE(obj$control$compute_viable_fitness)) {
    message("assembler> Computing viable bounds")
    obj$community <- community_viable_bounds(obj$community)
  }
  obj$done <- is.null(obj$community$bounds)
  obj <- assembler_append_history(obj)
  obj
}

assembler_births <- function(obj) {
  ## TODO: Change to:
  ## assembler_new_types <- function(obj) {
  ##   control <- obj$control
  ##   if (control$birth_type == "maximium") {
  ##     assembler_new_types_maximum_fitness(obj$community, control)
  ##   } else if (control$birth_type == "stochastic") {
  ##     assembler_new_types_stochastic(obj$community, control)
  ##   } else {
  ##     stop("Unknown birth type ", control$birth_type)
  ##   }
  ## }
  ## assembler_new_types_maximum_fitness <- function(sys, control) {
  ##   empty <- function(sys, m=NULL) {
  ##     ret <- trait_matrix(numeric(0), sys$trait_names)
  ##     ret <- copy_attributes(m, ret, exclude=c("dim", "dimnames"))
  ##     attr(ret, "done") <- TRUE
  ##     ret
  ##   }
  ##   eps_too_close <- control$eps_too_close
  ##   ...
  ## }
  to_add <- obj$births(obj$community)

  ## Now we add things:
  if (nrow(to_add) > 0) {
    message_new_types(to_add)
    i <- should_move(obj, to_add)
    if (length(i) == 1L) {
      obj$community <- assembler_births_try_move(obj$community, to_add, i)
    } else {
      rain <- initial_seed_rain(attr(to_add, "fitness"), obj$community)
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

assembler_births_try_move <- function(community, to_add, i) {
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
                    initial_seed_rain(w, community))
  } else {
    message("assembler[births]> Move was successful")
    community <- test
  }
  community
}

assembler_run_model <- function(obj) {
  run_type <- obj$control$run_type
  message(sprintf("assembler> Running model (%s)", run_type))
  message_community_state(obj$community)
  run <- switch(run_type,
                single=community_run,
                to_equilibrium=community_run_to_equilibrium,
                stop("Unknown run type ", run_type))
  obj$community <- run(obj$community)
  obj
}

assembler_deaths <- function(obj) {
  obj$community <- obj$deaths(obj$community)
  obj
}

assembler_restore <- function(obj, prev) {
  message("assembler> Restoring previous community and history")
  obj$history <- prev
  obj$community <- restore_community(last(obj$history), recompute=TRUE)
  obj$done <- is.null(obj$community$bounds)
  obj
}

assembler_prepare_fitness <- function(obj) {
  if (obj$control$birth_type == "maximum") {
    message("assembler> Computing approximate fitness")
    obj$community <- community_prepare_fitness(obj$community)
    obj$community <- community_prepare_approximate_fitness(obj$community)
  }
  obj
}

assembler_append_history <- function(obj) {
  obj <- assembler_prepare_fitness(obj)
  obj$history <- c(obj$history, list(obj$community))
  filename <- obj$filename
  if (!is.null(filename)) {
    ok <- try(saveRDS(obj$history, filename))
    if (inherits(ok, "try-error")) {
      warning("History saving has failed",
              immediate.=TRUE, call.=FALSE)
    }
  }
  obj
}

assembler_step <- function(obj) {
  message(sprintf("assembler> *** Assembler: step %d, (%d strategies), %s",
                  length(obj$history), length(obj$community),
                  obj$run_type))
  message_community_state(obj$community)
  if (!obj$done) obj <- assembler_births(obj)
  if (!obj$done) obj <- assembler_run_model(obj)
  if (!obj$done) obj <- assembler_deaths(obj)
  if (!obj$done) obj <- assembler_append_history(obj)
  obj
}

##' @export
length.assembler <- function(x) {
  length(x$history)
}

should_move <- function(obj, to_add) {
  ret <- integer(0)
  if (nrow(to_add) == 1L &&
      length(obj$community) > 0L &&
      obj$control$run_type == "to_equilibrium") {
    ## NOTE: not rescaled by obj$bounds?
    i <- closest_log(to_add, obj$community$traits)
    if (attr(i, "distance") < obj$control$birth_move_tol) {
      ret <- i
    }
  }
  ret
}

## Helper function for printing community state.  Not fast!
message_community_state <- function(community, prefix="assembler> ",
                                    header="*** Traits: ") {
  str <- format_community_state(community$traits,
                                community$seed_rain,
                                prefix)
  message(paste0(prefix, header))
  message(paste(str, collapse="\n"))
}

message_new_types <- function(to_add) {
  header <- "*** Proposed new type(s):"
  prefix <- "assembler[births]> "
  rain <- exp(attr(to_add, "fitness"))
  str <- format_community_state(to_add, rain, prefix)
  message(paste0(prefix, header))
  message(paste(str, collapse="\n"))
}

format_community_state <- function(traits, seed_rain, prefix=NULL) {
  m <- cbind(as.data.frame(traits), seed_rain=prettyNum(seed_rain))
  if (nrow(m) == 0) {
    m[1,] <- rep("<empty>", ncol(m))
    rownames(m) <- ""
  }
  str <- capture.output(print(m))
  if (!is.null(prefix)) {
    str <- paste0(prefix, str)
  }
  str
}
