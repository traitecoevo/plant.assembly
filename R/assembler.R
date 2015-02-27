assembler <- function(community, control=NULL, filename=NULL, prev=NULL) {
  control <- assembler_control(control)

  ret <- list(history=list(),
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
                   vcv=NULL, # will be needed...
                   ## community_deaths:
                   check_inviable=TRUE,
                   dead_seed_rain=1e-3,
                   ## births_maximum_fitness:
                   eps_too_close=1e-3)
  if (identical(control[["birth_type"]], "stochastic")) {
    defaults$run_type <- "single"
  }

  control <- as.list(control)
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown control parameters ", paste(extra, collapse=", "))
  }
  ret <- modifyList(defaults, control)
  if (ret$birth_type == "maximum" && ret$run_type != "to_equilibrium") {
    stop("Must use 'to_equilibrium' run type with maximum fitness births")
  }
  if (ret$birth_type == "stochastic" && is.null(ret$vcv)) {
    stop("vcv must be provided")
  }
  ret
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

assembler_restore <- function(obj, prev) {
  message("assembler> Restoring previous community and history")
  obj$history <- prev
  obj$community <- restore_community(last(obj$history), recompute=TRUE)
  obj$done <- is.null(obj$community$bounds)
  obj
}

assembler_prepare_fitness <- function(obj) {
  message("assembler> Computing ode times")
  obj$community <- community_prepare_fitness(obj$community)
  if (obj$control$birth_type == "maximum") {
    message("assembler> Computing approximate fitness")
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
                  obj$control$run_type))
  message_community_state(obj$community)
  if (!obj$done) obj <- assembler_births(obj)
  if (!obj$done) obj <- assembler_run_model(obj)
  if (!obj$done) obj <- assembler_deaths(obj)
  if (!obj$done) obj <- assembler_append_history(obj)
  obj
}

assembler_set_traits <- function(obj, traits, seed_rain=NULL) {
  if (length(obj$community) > 0L) {
    stop("This is for an empty community only")
  }
  message("assembler> Setting traits")
  obj$community <- community_add(obj$community, traits, seed_rain)
  message_community_state(obj$community)
  obj <- assembler_run_model(obj)
  obj <- assembler_deaths(obj)
  obj <- assembler_append_history(obj)
  obj
}

assembler_run <- function(obj, nsteps) {
  for (i in seq_len(nsteps)) {
    if (obj$done) {
      message(sprintf("assembler> Assembly completed after %d steps", i))
      break
    }
    obj <- assembler_step(obj)
  }
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
