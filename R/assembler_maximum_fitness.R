make_births_maximum_fitness <- function(control) {
  ## This bit of weirdness exists so that extra information from the
  ## fitness search, such as approxiate fitness points, gets copied
  ## back so we can work with it.
  empty <- function(sys, m=NULL) {
    ret <- trait_matrix(numeric(0), sys$trait_names)
    ret <- copy_attributes(m, ret, exclude=c("dim", "dimnames"))
    attr(ret, "done") <- TRUE
    ret
  }
  eps_too_close <- control$eps_too_close

  function(sys) {
    if (is.null(sys$bounds)) {
      return(empty(sys))
    }

    eps_equilibrium <- sys$parameters$control$equilibrium_eps
    eps_positive_fitness <- eps_equilibrium
    ret <- find_max_fitness(sys, eps_too_close)

    if (attr(ret, "fitness") < eps_positive_fitness) {
      message("births[max]> Best point had nonpositive fitness: ",
              attr(ret, "fitness"))
      return(empty(sys, ret))
    }

    if (length(sys) > 0L) {
      i <- closest_log(drop(ret), sys$traits, sys$bounds)
      if (attr(i, "distance") < eps_too_close) {
        message("births[max]> Best point too close to existing")
        return(empty(sys, ret))
      }
    }

    ret
  }
}

## This is fundamentally a really hard problem because we want to
## chase a lot of data back to the main object.  In particular:
##
## * We might not have an up-to-date cohort schedule or ode times: we
##   want to set them.
## * We want to get the full approximate fitness landscape
find_max_fitness <- function(sys, eps_too_close=1e-3) {
  bounds <- check_bounds(sys$bounds)
  if (nrow(bounds) == 1L) {
    find_max_fitness_1d(sys, eps_too_close)
  } else {
    find_max_fitness_2d(sys, eps_too_close)
  }
}

## NOTE: In the 1d case we don't use the optimise approach because it
## might miss local peaks.  Instead we construct an approximate
## landscape and look around the highest point.
find_max_fitness_1d <- function(sys, eps_too_close) {
  bounds <- check_bounds(sys$bounds, finite=TRUE)

  ## This should be cheap, but will be lost if the calling function
  ## didn't arrange it:
  sys <- community_prepare_approximate_fitness(sys)
  fitness_approximate <- community_fitness_approximate(sys)

  xx <- seq_log_range(sys$bounds, 500)
  yy <- fitness_approximate(xx)
  i <- which.max(yy)

  ## If we want to polish this point a bit, we could optimise over the
  ## actual fitness function or the approximate; for now I'm using the
  ## approximate as this will be much faster and the optimum should
  ## actually lie in that range.
  r <- xx[c(max(1, i - 1), min(i + 1, length(xx)))]
  opt <- optimise(fitness_approximate, r, maximum=TRUE)

  ret <- trait_matrix(opt$maximum, sys$trait_names)
  attr(ret, "fitness") <- opt$objective
  ret
}

find_max_fitness_2d <- function(sys, eps_too_close, tol=1e-2) {
  stop("Not yet implemented")
}
