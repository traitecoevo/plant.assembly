##' Initialise a community assembly object
##'
##' Control options affect how the assembly proceeds.
##'
##' @title start community assembly process.
##' @param community A \code{community} object.
##' @param control A control object affecting how the assembly proceeds,
##' see \code{assembler_control} for details,
##' @param filename Location to save results.
##' @param prev History from previously run assembly, used to restart.
##' @return An \code{assembler} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
assembler_start <- function(community, control=assembler_control(control), filename=NULL, prev=NULL) {
  
  ret <- list(history=list(),
              done=FALSE,
              control=control,
              filename=filename)

  if (!is.null(filename) && file.exists(filename)) {
    if (is.null(prev)) {
      prev <- readRDS(filename)
    } else if (!isTRUE(all.equal(readRDS(filename), prev))) {
      stop("filename present and prev given, but don't agree!")
    }
  }

  if (is.null(prev)) {
    ret <- assembler_initialise(ret, community)
  } else {
    ret <- assembler_restore(ret, community, prev)
  }

  class(ret) <- "assembler"
  ret
}

##' Controls how community assembly works.
##'
##' Returns a list. Passing in a list of value via \code{
##' control} will override the defaults. Options include
##' run_type determines whether population is stepped to
##' demographic equilibrium ("to_equilibrium") or not ("single").
##' "birth_type" determines sampling of new types -- "stochastic" or
##' "maximum" (on fitness peak). With "stochastic" births,
##' "n_mutants" and "n_immigrants" determine the frequency of
##' resident mutations and immigrations from global pool.
##' "vcv" is variance-covariance matrix for mutations.
##' If "birth_move_tol" is trait distance
##' within which we attempt to move an existing resident rather introduce
##' a new type (this helps reduce the number of types).
##' "compute_viable_fitness" asks whether to check bounds of viable
##' trait space. "check_positive" determines whether the fitness of an
##' invader is checked before it is introduced. If
##' "check_inviable" causes dead residents to be removed when birth rate
##' drops below "dead_birth_rate".
##' "eps_too_close" is tolerance in trait values when searching for maxima.
##'
##' @title Options controllings community assembly process.
##' @param control A list of values to modify from defaults.
##' @return A list with elements run_type, birth_type, birth_move_tol,
##' compute_viable_fitness, n_mutants, n_immigrants, check_positive,
##' vcv, check_inviable, dead_birth_rate, eps_too_close
##' @author Rich FitzJohn, Daniel Falster
##' @export
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
                   dead_birth_rate=1e-3,
                   ## births_maximum_fitness:
                   eps_too_close=1e-3,
                   min_birth_rate_initial = 1e-3,
                   max_birth_rate_initial = 500,
                  ## This magic number is not great, but needed to 
                  ## prevent suggesting adding 1e8 as the birth 
                  ## rate (can happen!).  This should be "quite 
                  ## large", but obviously that number depends on
                  ## the situation.
                  equilibrium_eps = 0.001
                  )
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
  plant_log_assembler("Starting empty assembler")
  obj$community <- community
  if (isTRUE(obj$control$compute_viable_fitness)) {
    plant_log_assembler("Computing viable bounds")
    obj$community <- community_viable_bounds(obj$community)
  }
  obj$done <- is.null(obj$community$bounds)
  obj <- assembler_append_history(obj)
  obj
}

assembler_run_model <- function(obj) {
  run_type <- obj$control$run_type
  plant_log_assembler(sprintf("Running model (%s)", run_type))
  plant_log_assembler_state(obj$community)
  run <- switch(run_type,
                single=community_run,
                to_equilibrium=community_run_to_equilibrium,
                stop("Unknown run type ", run_type))

  obj$community <- run(obj$community)
  obj
}

assembler_restore <- function(obj, community, prev) {
  plant_log_assembler("Restoring previous community and history")

  if (length(community) > 0L) {
    stop("community must be an empty placeholder")
  }
  ## Check that this looks legit:
  x <- last(prev$history)
  if (!isTRUE(all.equal(community$parameters, x$parameters))) {
    stop("Restoring community: parameters don't match")
  }
  if (!all(x$bounds[, "lower"] >= community$bounds[, "lower"] &
           x$bounds[, "upper"] <= community$bounds[, "upper"])) {
    stop("Restoring community: incompatible bounds")
  }
  plant_log_assembler(sprintf("...restored %d steps of history",
                              length(prev$history)))
  prev
}

assembler_prepare_fitness <- function(obj) {
  plant_log_assembler("Computing ode times")
  obj$community <- community_prepare_fitness(obj$community)
  if (obj$control$birth_type == "maximum") {
    plant_log_assembler("Computing approximate fitness")
    obj$community <- community_prepare_approximate_fitness(obj$community)
  }
  obj
}

assembler_append_history <- function(obj) {
  obj <- assembler_prepare_fitness(obj)
  obj$history <- c(obj$history, list(obj$community))
  if (isTRUE(obj$done)) {
    attr(obj$history, "done") <- TRUE
  }
  filename <- obj$filename
  if (!is.null(filename)) {
    ok <- try(saveRDS(obj, filename))
    if (inherits(ok, "try-error")) {
      warning("History saving has failed",
              immediate.=TRUE, call.=FALSE)
    }
  }
  obj
}

assembler_step <- function(obj) {
  plant_log_assembler(sprintf("step %d, (%d strategies), %s",
                  length(obj$history), length(obj$community),
                  obj$control$run_type))
  plant_log_assembler_state(obj$community)
  if (!obj$done) obj <- assembler_births(obj)
  if (!obj$done) obj <- assembler_run_model(obj)
  if (!obj$done) obj <- assembler_deaths(obj)
  assembler_append_history(obj)
}

assembler_set_traits <- function(obj, traits, birth_rate=NULL) {
  if (length(obj$community) > 0L) {
    stop("This is for an empty community only")
  }
  plant_log_assembler("Setting traits")
  obj$community <- community_add(obj$community, traits, birth_rate)
  plant_log_assembler_state(obj$community)
  obj <- assembler_run_model(obj)
  obj <- assembler_deaths(obj)
  obj <- assembler_append_history(obj)
  obj
}

##' Take specified number of steps with assembler.
##'
##'
##' @title Take specified number of steps with assembler.
##' @param obj An  \code{assembler} object.
##' @param nsteps Number of steps
##' @return An \code{assembler} object.
##' @author Rich FitzJohn, Daniel Falster
##' @export
assembler_run <- function(obj, nsteps) {
  for (i in seq_len(nsteps)) {
    if (obj$done) {
      plant_log_assembler(sprintf("Assembly completed after %d steps",
                                  length(obj$history)))
      break
    }
    obj <- assembler_step(obj)
  }
  obj
}

##' Returns number taken by an assembler object
##'
##'
##' @title Returns number taken by an \code{assembler} object
##' @param x An  \code{assembler} object.
##' @return Length of history in assembler object.
##' @author Rich FitzJohn, Daniel Falster
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
##
## TODO: This might need a trailing newline added after?  I'm getting
## no trailing newline.  That's a bit surprising and possibly a bug in
## loggr.
plant_log_assembler_state <- function(community) {
  str <- format_community_state(community$traits,
                                community$birth_rate,
                                NULL)
  msg <- paste(c("*** Traits:", str), collapse="\n")
  plant_log_assembler(msg,
                      traits=community$traits,
                      birth_rate=community$birth_rate)
}

format_community_state <- function(traits, birth_rate, prefix=NULL) {
  m <- cbind(as.data.frame(traits, stringsAsFactors=FALSE),
             birth_rate=I(prettyNum(birth_rate)))
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
