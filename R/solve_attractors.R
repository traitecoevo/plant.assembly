## Functions *related* to assembly, but not actually doing it.

##' Compute max growth rate of a given set of values for a trait.
##' This is the log of per-capita seed production (i.e., fitness).
##'
##' Only works in one dimension
##' @title Compute Max Growth Rate
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param values Values to compute maximum growth rate for
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param schedule \code{CohortSchedule} to use, or \code{NULL} to
##' generate a hopefully reasonable schedule.
##' @author Rich FitzJohn
##' @export
max_growth_rate <- function(trait, values, p, schedule=NULL) {
  fitness_landscape_empty(trait, values, p, schedule)
}

##' Compute the carrying capacity (equilibrium per-capita seed
##' production) for a set of values of a trait.  Each is considered in
##' isolation.
##'
##' @title Carrying Capacity
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param values Values to compute maximum growth rate for
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param birth_rate Initial seed rain (optional)
##' @param parallel Use multiple processors?
##' @author Rich FitzJohn
##' @export
carrying_capacity <- function(trait, values, p, birth_rate=1,
                              parallel=FALSE) {
  f <- function(x) {
    carrying_capacity1(trait, x, p, birth_rate)
  }
  if (parallel) {
    unlist(parallel::mclapply(values, f, mc.preschedule=FALSE))
  } else {
    unlist(lapply(values, f))
  }
}

carrying_capacity1 <- function(trait, value, p, birth_rate = 1) {
  p <- p$copy()
  p$clear()
  ## Use
  p <- expand_parameters(trait, value, p)
  p$birth_rate <- birth_rate
  if (is.null(schedule)) {
    schedule <- default_cohort_schedule(p)
  }
  res <- equilibrium_birth_rate(p)
  warning("Please fix seed rain mean", immediate. = TRUE)
  rowMeans(res$birth_rate)
}


##' Compute region of positive fitness.  This will be the values where
##' fitness is approximately zero.
##'
##' @title Compute Region of Positive Fitnes
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param bounds 2D vector specifing range within which to search
##' @param value Initial value - must have positive fitness itself!
##' If not given, then the value from the default strategy within
##' \code{p} is used.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @param dx Amount to step the trait.  If \code{log_scale} is
##' \code{TRUE}, this is on a log scale.
##' @param find_max_if_negative If the starting value has negative
##' fitness, should we search for the maximum value?
##' @export
##' @author Rich FitzJohn

positive <- function(f, x, dx, lower=-Inf, upper=Inf, eps=1e-3) {
  b <- plant:::positive_1d_bracket(f, x, dx, lower, upper)

  # Find lower root. If no root exists within that range, take
  # lower bound
  if(prod(b$lower$fx[1:2]) < 0){
    ## The suppressWarnings here is about conversion from -Inf to a
    ## very small number.
    lower <- suppressWarnings(uniroot(f, b$lower$x,
                                      f.lower=b$lower$fx[[1]],
                                      f.upper=b$lower$fx[[2]],
                                      tol=eps)$root)
  } else {
    lower <- b$lower$x[[1]]
  }
  # Find upper root. If no root exists within that range, take
  # upper bound
  if(prod(b$upper$fx[1:2]) < 0){
    upper <- uniroot(f, b$upper$x,
                   f.lower=b$upper$fx[[1]], f.upper=b$upper$fx[[2]],
                   tol=eps)$root
  } else {
    upper <- b$upper$x[[2]]
  }

  c(lower, upper)
}

## This is a multidimensional version of positive.  It's a hack for
## now.
##
## Unlike positive, it *requires* bounds.
## It takes n_start, which is the *per dimension* number of points to
## start with.  n_total is the total number of points to test.  It
## also takes a starting point.  The idea is that starting point has
## positive fitness so we can use it as a root-creating device.
positive2 <- function(f, x, lower, upper, n_total=200) {
  delaunay_init()
  is_nonnegative <- function(y) y >= 0.0
  pts <- lapply(seq_along(lower), function(i) c(lower[[i]], upper[[i]]))
  m0 <- rbind(x, unname(as.matrix(do.call("expand.grid", pts))))
  y0 <- f(m0)
  delaunay_run_map(m0, f, is_nonnegative, values=y0,
                   n_total=n_total, exploit=50)
}

#' Solve for equilbrium community with given values of specified trait
#'
#' Find seed rain and cohort schedule for equilbrium community with
#' given traits. This is point at which each resident seed returns
#' on average a single seed.
#' @param trait name of trait
#' @param values vector of trait values for community
#' @param p Parameters object to use.  Importantly, the
#' \code{strategy_default} element gets used here.
#' @param birth_rate vector of initial seed rains for community
#' @author Daniel Falster
#' @export
get_equilibrium_community <- function(trait, values, p, birth_rate=NULL) {
  p <- p$copy()
  p$clear()
  p <- expand_parameters(trait, values, p, FALSE)
  if (!is.null(birth_rate)) {
    if (length(values) != length(birth_rate)) {
      stop("incorrect vector lengths")
    }
    p$birth_rate <- birth_rate
  }
  res <- equilibrium_birth_rate2(p)
  ## Take the *final*, not the mean value: this is important for
  ## assembly.  TODO: make this change elsewhere too.
  p$birth_rate <- unname(res$birth_rate[,"out"])
  list(p=p, schedule=res$schedule)
}

#' Find evolutionary attractor for single species and trait.
#'
#' Find evolutionary attractor for single species and trait.
#' This is point at which selection gradient equals zero.
#' Currently solved using \code{uniroot}
#' @param trait name of trait
#' @param p Parameters object to use.  Importantly, the
#' \code{strategy_default} element gets used here.
#' @param bounds a vector containing the end-points of the
#' interval to be searched for the root.
#' @param ... set verbpse=TRUE for verbose output, birth_rate=value
#' gives starting values when solving for demographic equilibrium
#' @param tol the desired accuracy (convergence tolerance).
#' @param edge_ok Is it (not) an error if we end up on the edge of the
#' viable bounds?
#' @author Daniel Falster
#' @export
#' @return A species object, with trait, seed rain and cohort schedule
#' information.
community_solve_singularity_1D <- function(community, bounds = NULL, tol = 1e-04, ...,
                                edge_ok = TRUE) {

  plant_log_assembler(
    sprintf("Solving 1D attractor for %s", community$trait_names))
                                  
  f <- function(x) {
    out <- 
      community %>%
      community_add(plant::trait_matrix(x, "lma")) %>%
      community_run_to_equilibrium() %>%
      community_selection_gradient_1d()
    
    ret <- out$selection_gradient
    
    # Add extra details so we can access these later
    attr(ret, "community") <- out

    ret
  }
    
  if(is.null(bounds)) 
    bounds <- community$bounds

  lower <- bounds[[1]]
  upper <- bounds[[2]]

  f_lower <- f(lower)
  f_upper <- f(upper)

  ## This is the exact condition used by uniroot:
  failed <- !isTRUE(as.vector(sign(f_lower) * sign(f_upper) <= 0))

  if (failed) {
    msg <- paste("Bounds do not include attractor: taking",
                 if (f_lower < 0) "lower" else "upper")
    if (edge_ok) {
      warning(msg)
    } else {
      stop(msg)
    }

    if (f_lower < 0) {
      res <- list(root = lower, f.root = f_lower)
    } else {
      res <- list(root = upper, f.root = f_upper)
    }
  } else {
    res <- uniroot(f,
      lower = lower, upper = upper,
      f.lower = f_lower, f.upper = f_upper, tol = tol
      )
  }

  # We're actualy going to 
  community_out <- attr(res$f.root, "community")
  if(community_out$traits != res$root) {
    stop("community_solve_singularity_1D: Hmm, these values should be equla. Better quit now.")
  }

  plant_log_assembler(
    sprintf("Solved! 1D attractor for %s is %s", community$trait_names, as.numeric(res$root))
  )

  community_out
}

#' Calculates selection gradient in a single-species community
#' with given trait value
#'
#' Adds selection gradient to a single-species community
#' with given trait value. This is derivative of fitness
#' with respect to trait value. You should first solv for
#' using \code{community_run_to_equilibrium}
#' @param community community object to use. 
#' @param dx Interval over which derivative is calculated
#' @param log_scale (currently disabled) Determines whether derivative is taken
#' with respect to raw or log-transformed x values. The latter
#' is useful when x is log-normally distributed.
#' @author Daniel Falster
#' @export
#' @return a community with selection gradient added.
community_selection_gradient_1d <- function(community, dx=1e-04,
                                log_scale=TRUE) {

  msg <- sprintf("Calculating selection gradient for %s = %f",community$trait_names, community$traits)
  plant_log_assembler(msg)
  
  trait_names <- community$trait_names
  # get points needed for gradient
  points <- community$traits[,trait_names] * c(1, 1 + dx)
  # format points
  traits <- plant::trait_matrix(points, trait_names)
  # calculate fitness
  ff <- community_fitness(community, traits)  
  # caluclate gradient using forward difference
  ret <- (ff[2] - ff[1]) / dx

  # alternative method using grader
  f  <- function(x) {
    traits <- plant::trait_matrix(x, "lma")
    community_fitness(community, traits)
  } 
  #grader::gradient(f, community$traits[, trait])

  msg <- sprintf("Solved! Selection gradient for %s = %f is %s", community$trait_names, community$traits, ret)
  plant_log_assembler(msg)

  community[["fitness"]] <- ff[1]
  community[["selection_gradient"]] <- ret
  community
}

##' Find point of maximum fitness in empty fitness landscape within a
##' specified range.
##'
##' @title Find point of maximum fitness within some range.
##' @param trait Name of the trait (e.g., \code{"lma"})
##' @param bounds Two element vector specifing range within which to
##' search
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @export
##' @author Daniel Falster, Rich FitzJohn
max_fitness <- function(trait, p, bounds=NULL, log_scale=TRUE) {
  if(length(trait) > 1) {
    stop("Doesn't yet support multiple traits")
  }
  ## These are unlikely to be very good in general.  Bounds must be
  ## finite.
  if (is.null(bounds)) {
    bounds <- c(1e-5, 1e3)
  }
  if (log_scale) {
    bounds <- log(bounds)
    f <- function(x) max_growth_rate(trait, exp(x), p)
  } else {
    f <- function(x) max_growth_rate(trait, x, p)
  }

  out <- suppressWarnings(optimise(f, interval=bounds, maximum=TRUE, tol=1e-3))
  structure(if (log_scale) exp(out$maximum) else out$maximum,
            fitness=out$objective)
}
