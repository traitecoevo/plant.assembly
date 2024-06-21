plant_log_eq <- function(...) {
  plant_log_info(..., routine = "equilibrium")
}


##' Run system to offspring arrival equilibrium
##'
##' @title Run system to offspring arrival equilibrium
##' @param p A \code{Parameters} object
##' @param ctrl Control object
##' @return A Parameters object, with offspring arrival and node schedule
##' elements set.
##' @export
##' @author Rich FitzJohn
equilibrium_birth_rate <- function(p, ctrl) {
  solver <- ctrl$equilibrium_solver_name
  plant_log_info(sprintf("Solving offspring arrival using %s", solver),
                 routine = "equilibrium", stage = "start", solver = solver)
  switch(solver,
         iteration = equilibrium_birth_rate_iteration(p, ctrl = ctrl),
         hybrid = equilibrium_birth_rate_hybrid(p, ctrl = ctrl),
         nleqslv = equilibrium_birth_rate_solve(p, ctrl = ctrl, solver = "nleqslv"),
         dfsane = equilibrium_birth_rate_solve(p, ctrl = ctrl, solver = "dfsane"),
         stop("Unknown solver ", solver))
}

## This is the simplest solver: it simply iterates the outgoing offspring
## produced as incoming offspring arrival.  No attempt at projection is made.
equilibrium_birth_rate_iteration <- function(p, ctrl) {
  
  check <- function(x_in, x_out, eps, verbose) {
    achange <- x_out - x_in
    rchange <- 1 - x_out / x_in
    ## eps > 0 && # <- this was a precondition - seems odd.
    converged <- all(abs(achange) < eps | abs(rchange) < eps)
    if (converged) {
      achange <- max(abs(achange))
      rchange <- max(abs(rchange))
      fmt <- "Reached target accuracy (delta %2.5e, %2.5e < %2.5e eps)"
      plant_log_eq(sprintf(fmt, achange, rchange, eps),
                   stage="converged", achange=achange, rchange=rchange)
    }
    converged
  }

  eps <- ctrl$equilibrium_eps
  verbose <- ctrl$equilibrium_verbose
  
  birth_rates <- purrr::map_dbl(p$strategies, ~ purrr::pluck(., "birth_rate_y"))

  runner <- make_equilibrium_runner(p, ctrl = ctrl)
  
  for (i in seq_len(ctrl$equilibrium_nsteps)) {
    offspring_production <- runner(birth_rates)
    converged <- check(birth_rates, offspring_production, eps, verbose)
    birth_rates <- offspring_production
    if (converged) {
      break
    }
  }

  # TODO: revisit 'gross' behaviour in cleanup utility
  equilibrium_runner_cleanup(runner, converged)
}

equilibrium_birth_rate_solve <- function(p, ctrl = scm_base_control(),
                                         solver="nleqslv") {
  try_keep <- ctrl$equilibrium_solver_try_keep
  logN <- ctrl$equilibrium_solver_logN
  min_offspring_arriving <- 1e-10 # TODO: should also be in the controls?

  plant_log_eq(paste("Solving offspring arrival using", solver),
               stage="start", solver=solver)

  birth_rates <- purrr::map_dbl(p$strategies, ~ purrr::pluck(., "birth_rate_y"))
  runner <- make_equilibrium_runner(p,ctrl =ctrl)

  ## First, we exclude species that have offspring arrivals below some minimum
  ## level.
  to_drop <- birth_rates < min_offspring_arriving
  if (any(to_drop)) {
    i_keep <- which(!to_drop)
    msg <- sprintf("Species %s extinct: excluding from search",
                   paste(which(to_drop), collapse=" & "))
    plant_log_eq(msg, stage="drop species", drop=which(to_drop))
    offspring_arriving_full <- birth_rates
    birth_rates <- offspring_arriving_full[i_keep]

    runner_full <- runner
    runner <- function(x) {
      x_full <- rep(0, length(offspring_arriving_full))
      x_full[i_keep] <- x
      runner_full(x_full)[i_keep]
    }
  }

  ## Then see if any species should be retained:
  if (try_keep) {
    ans <- runner(birth_rates)
    keep <- unname(ans >= birth_rates)

    msg <- sprintf("Keeping species %s",
                   paste(which(!to_drop)[keep], collapse=", "))
    plant_log_eq(msg, stage="keep species", keep=which(!to_drop)[keep])
  } else {
    keep <- rep(FALSE, length(p$strategies))
  }

  ## TODO: This is annoying, but turns out to be a problem for getting
  ## the solution working nicely.
  max_offspring_arriving <- pmax(birth_rates * 100, 10000)
  target <- equilibrium_birth_rate_solve_target(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving,
                                               ctrl$equilibrium_verbose)
  x0 <- if (logN) log(birth_rates) else birth_rates

  tol <- ctrl$equilibrium_eps
  ## NOTE: Hard coded minimum of 100 steps here.
  maxit <- max(100,
               p$control$equilibrium_nsteps)
  sol <- util_nlsolve(x0, target, tol = tol, maxit = maxit, solver = solver)
  res <- equilibrium_runner_cleanup(runner, attr(sol, "converged"))
  attr(res, "sol") <- sol
  res
}

## The idea is to use rounds of iteration to try and push the
## system into the basin of attraction of the stable equilibrium.  The
## final approach is slow so use a root-finding approach there.
## However, if we are *not* in the basin of attraction the root finder
## will happily select zero offspring arrivals for species that are not
## zeros.  So after running a round with the solver, check any species
## that were zeroed to make sure they're really dead.
equilibrium_birth_rate_hybrid <- function(p, ctrl) {
  attempts <- ctrl$equilibrium_nattempts

  ## Then expand this out so that we can try alternating solvers
  solver <- rep(c("nleqslv", "dfsane"), length.out=attempts)

  for (i in seq_len(attempts)) {
    eq_solution_iteration <- equilibrium_birth_rate_iteration(p, ctrl = ctrl)
    
    converged_it <- isTRUE(attr(eq_solution_iteration, "converged"))
    msg <- sprintf("Iteration %d %s",
                   i, if (converged_it) "converged" else "did not converge")
    plant_log_eq(msg, step="iteration", converged = converged_it, iteration=i)

    eq_solution <- try(
      equilibrium_birth_rate_solve(
        eq_solution_iteration, 
        ctrl = ctrl, 
        solver = solver[[i]]
        )
      )

    converged_sol <- isTRUE(attr(eq_solution, "converged"))

    msg <- sprintf("Solve %d %s",
                    i, if (converged_sol) "converged" else "did not converge")
    plant_log_eq(msg, step="solve", converged=converged_sol, iteration=i)

    if (converged_sol) {

      # check species with zero eq. birth rate are truly unviable.
      extinct = purrr::map_lgl(eq_solution$strategies, function(s) s$birth_rate_y == 0.0)

      if (any(extinct)) {
        plant_log_eq("Checking species driven to extinction")
        
        ## Add extinct species back at extremely low density and make sure
        ## that this looks like a legit extinction.
        p_check <- eq_solution
        p_check$strategies[extinct]$birth_rate_y <- ctrl$equilibrium_extinct_birth_rate

        res <- run_scm(p_check)

        # `next` breaks the loop iterating over solutions and does not return `eq_solution`
        if (any(res$offspring_production[extinct] > ctrl$equilibrium_extinct_birth_rate)) {
          plant_log_eq("Solver drove viable species extinct: rejecting")
          next
        }
      }
      plant_log_eq("Accepting solution via solver")
      return(eq_solution)
    }
  }

  ## This one should be a warning?
  plant_log_eq("Repeated rounds failed to find optimum; returning solution from equilibrium_birth_rate_iteration")
  return(eq_solution_iteration)
}

## Another layer of runner for the solver code:
equilibrium_birth_rate_solve_target <- function(runner, keep, logN,
                                               min_offspring_arriving, max_offspring_arriving,
                                               verbose) {
  force(runner)
  force(keep)
  force(logN)
  force(min_offspring_arriving)
  force(max_offspring_arriving)
  function(x, ...) {
    if (logN) {
      x <- exp(x)
    }
    ## TODO: most of the plant_log_eq things here should be DEBUG not INFO?
    ## Avoid negative offspring arrivals:
    x[x < min_offspring_arriving & keep] <- min_offspring_arriving
    x[x < min_offspring_arriving & !keep] <- 0.0
    if (!any(x > 0)) {
      plant_log_eq("All species extinct?")
    }
    too_high <- x > max_offspring_arriving
    if (any(too_high)) {
      plant_log_eq(sprintf("Truncating offspring arrival of species %s",
                           paste(which(too_high), collapse=", ")))
      if (length(max_offspring_arriving) == 1L) {
        max_offspring_arriving <- rep(max_offspring_arriving, length.out=length(x))
      }
      x[too_high] <- max_offspring_arriving[too_high]
    }
    xout <- unname(runner(x))

    xout[keep] <- xout[keep] / x[keep] - 1.0
    i <- !keep & x > 0
    xout[i] <- xout[i] - x[i]

    ## !keep & x <= 0 will now  be zero
    if (any(xout[!keep & x <= 0] != 0.0)) {
      warning("This is not expected", immediate.=TRUE)
    }
    xout
  }
}

## Support code:
make_equilibrium_runner <- function(p, ctrl) {
  pretty_num_collapse <- function(x, collapse = ", ") {
    paste0("{", paste(prettyNum(x), collapse = collapse), "}")
  }

#  p <- validate(p)

  # default is about 10 ind.m-2
  large_offspring_arriving_change <- ctrl$equilibrium_large_birth_rate_change

  # traverse list of strategies and pull birth_rates
  last_arrival_rates <- purrr::map_dbl(
    p$strategies,
    ~ purrr::pluck(., "birth_rate_y")
  )

  default_schedule_times <- rep(
    list(p$node_schedule_times_default),
    length(last_arrival_rates)
  )

  last_schedule_times <- p$node_schedule_times
  history <- NULL
  counter <- 1L

  function(birth_rates) {
    # if a runner has diverged significantly in the last iteration then
    # reset the schedule to it's defaults and rebuild
    if (any(abs(birth_rates - last_arrival_rates) > large_offspring_arriving_change)) {
      p$node_schedule_times <- default_schedule_times
    }

    # set birth rates
    for (i in seq_along(p$strategies)) {
      p$strategies[[i]]$birth_rate_y <- birth_rates[i]
    }

    # update schedule - starts from prev. schedule so should be fast for fine scale resolution
    p_new <- build_schedule(p, ctrl = ctrl)

    # TODO: change behaviour of `build_schedule` to objects rather than attributes
    offspring_production <- attr(p_new, "offspring_production", exact = TRUE)


    # (ANDREW) Double arrow modifies counter in parent environment..
    # I don't love this approach to counting, as it introduces a side effect
    # to an otherwise functional programming design

    ## These all write up to the containing environment:
    p <<- p_new
    last_schedule_times <<- p_new$node_schedule_times
    history <<- c(history, list(c("in" = birth_rates, "out" = offspring_production)))

    msg <- sprintf(
      "eq> %d: %s -> %s (delta = %s)", counter,
      pretty_num_collapse(birth_rates),
      pretty_num_collapse(offspring_production),
      pretty_num_collapse(offspring_production - birth_rates)
    )
    plant_log_eq(msg,
      stage = "runner",
      iteration = counter,
      birth_rate = birth_rates,
      offspring_production = offspring_production
    )


    counter <<- counter + 1L

    # TODO: check why returning node schedule as an attribute
    # attr(offspring_production, "schedule_times") <- last_schedule_times

    offspring_production
  }
}

equilibrium_runner_cleanup <- function(runner, converged = TRUE) {
  # this pulls the history of the runner from the parent environment
  e <- environment(runner)
  if (is.function(e$runner_full)) {
    runner <- e$runner_full
    e <- environment(runner)
  }

  # the final solution has a recently built integration schedule
  p <- e$p

  # ANDREW:
  # I'm pretty sure this `p` has already been through `build_schedule`
  # but overloading the schedule times because that's what this used to do.
  p$node_schedule_times <- e$last_schedule_times

  attr(p, "progress") <- util_rbind_list(e$history)
  attr(p, "converged") <- converged
  p
}
