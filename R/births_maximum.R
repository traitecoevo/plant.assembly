community_new_types_maximum_fitness <- function(sys, control) {
  
  ## This bit of weirdness exists so that extra information from the
  ## fitness search, such as approxiate fitness points, gets copied
  ## back so we can work with it.
  empty <- function(sys, m=NULL) {
    ret <- trait_matrix(numeric(0), sys$trait_names)
    ret <- copy_attributes(m, ret, exclude=c("dim", "dimnames"))
    attr(ret, "done") <- TRUE
    ret
  }

  if (is.null(sys$bounds)) {
    ## TODO : can we delete `empty` function above
    browser()
    return(empty(sys))
  }

  ret <- find_max_fitness(sys, control)

  if (attr(ret, "fitness") < control$eps_fitness_invasion) {
    plant_log_max_fitness(paste0("Best point had nonpositive fitness: ",
                                 attr(ret, "fitness")))
    return(empty(sys, ret))
  }

  if (length(sys) > 0L) {
    i <- closest_log(drop(ret), sys$traits, sys$bounds)
    if (attr(i, "distance") < control$eps_too_close) {
      plant_log_max_fitness("Best point too close to existing")
      return(empty(sys, ret))
    }
  }

  ret
}

## This is fundamentally a really hard problem because we want to
## chase a lot of data back to the main object.  In particular:
##
## * We might not have an up-to-date cohort schedule or ode times: we
##   want to set them.
## * We want to get the full approximate fitness landscape
find_max_fitness <- function(sys, control) {

  plant_log_assembler("Finding maximum in fitness landscape")

  bounds <- check_bounds(sys$bounds)
  if (nrow(bounds) == 1L) {
    find_max_fitness_1D(sys, control$eps_too_close)
  } else {
    find_max_fitness_2d(sys, control$eps_too_close)
  }
}

## Simple function - takes existing maximum and then runs optim. Risks missing maximum
find_max_fitness_1D <- function(sys, eps_too_close) {

  # option 1 - use existing points
  i <- which.max(sys$fitness_points$fitness)
  xx <- sys$fitness_points[, 1, drop = TRUE]

  # todo: option 2 - fit 1D GP and use this to find max
  # xx <- seq_log_range(sys$bounds, 500)
  # yy <- fitness_approximate(xx)

  # Polish root by optnmising with actual fitness function
  ## first find range over which to look 
  r <- xx[c(max(1, i - 1), min(i + 1, length(xx)))]
  f <- sys$fitness_function
  
  opt <- optimise(f, r, maximum = TRUE)

  ret <- trait_matrix(opt$maximum, sys$trait_names)
  attr(ret, "fitness") <- opt$objective
  ret
}


find_max_fitness_2d <- function(sys, eps_too_close, tol=1e-2) {
  do_fit <- function(p) {
    f <- function(x) {
      community_fitness(sys, trait_matrix(x, sys$trait_names))
    }
    fit <- maximize_logspace(f, p, sys$bounds, tol)
    ret <- trait_matrix(fit$par, sys$trait_names)
    attr(ret, "fitness") <- fit$value
    ret
  }
  check <- function(fit, X) {
    j <- closest(log(drop(fit)), log(X), log(sys$bounds))
    w <- attr(fit, "fitness")
    d <- attr(j, "distance")
    plant_log_max_fitness(sprintf("\t...fitness: %s, distance: %s from %d",
                                  prettyNum(w), prettyNum(d), j))
    d > eps_too_close && w > 0.0
  }

  if (length(sys) == 0L) {
    p0 <- exp(rowMeans(log(sys$bounds)))
    ret <- do_fit(p0)
  } else {
    ret <- NULL
    X <- sys$traits
    gr_norm <- vnapply(sys$fitness_slopes, function(x) norm2(x$gr))
    idx <- order(gr_norm, decreasing=TRUE)
    attempts <- list()

    for (i in idx) {
      plant_log_max_fitness(paste0("Searching from species ", i))
      fit <- do_fit(X[i, ])
      attempts <- c(attempts, list(fit))
      if (check(fit, X)) {
        ret <- fit
        break
      }
    }
    if (is.null(ret)) {
      ## For want of a better thing to try:
      plant_log_max_fitness("Searching from the middle of occupied space")
      p0 <- exp(rowMeans(log(sys$bounds)))
      ret <- do_fit(p0)
      ## Note that we don't check this point: we'll do that in the
      ## make_births_maximum_fitness.
    }

    attr(ret, "attempts") <- ret
  }
  ret
}
