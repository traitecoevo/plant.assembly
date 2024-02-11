
#' Find point of maximum fitness in empty fitness landscape within a
##' specified range.
##'
##' @title Find point of maximum fitness within some range.
##'
##' @param bounds Two element vector specifying range within which to
##' search
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.
##' @param params Parameters for fundamental_fitness function
##' @param tol Tolerance used in the optimisation
##'
##' @importFrom stats optimise optim 
##' @export
##' @author Daniel Falster, Rich FitzJohn


solve_max_fitness <- function(bounds, params, log_scale = TRUE, tol = 1e-3){

  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)

  if (log_scale) {
    bounds[bounds[,1] == -Inf, 1] <- 0
    bounds <- log(bounds)

    ff <- exp
  } else {
    ff <- I
  }

  f <- function(x) fundamental_fitness(ff(trait_matrix(x, traits)), params)

  ret <- solve_max_worker(bounds, f, tol = 1e-3, outcome ="fitness")
  #browser()
  if (log_scale) {
    ret <- exp(ret)
  }

  return(ret)
}

solve_max_worker <- function(bounds, f, tol=1e-3, outcome) {
  if (length(rownames(bounds)) == 1L) {
    if (!all(is.finite(bounds))) {
      stop("Starting value did not have finite fitness; finite bounds required")
    }
    ## The suppressWarnings here is for warnings like:
    ##
    ## Warning message:
    ## In optimise(f, interval = bounds, maximum = TRUE, tol = tol) :
    ##   NA/Inf replaced by maximum positive value
    ##
    ## which is probably the desired behaviour here.
    out <- suppressWarnings(optimise(f, interval=bounds, maximum=TRUE, tol=tol))
    # browser()
    ret <- out$maximum
    attr(ret, outcome) <- out$objective

  } else {
    ## This is not very well tested, and the tolerance is not useful:
    out <- optim(rowMeans(bounds), f, method="L-BFGS-B",
                 lower=bounds[, "lower"], upper=bounds[, "upper"],
                 control=list(fnscale=-1, factr=1e10))

    ret <- out$value
    attr(ret, outcome) <- out$par
  }
  return(ret)
}


max_fitness <- function(bounds, p, log_scale=TRUE, tol=1e-3) {
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)

  if (log_scale) {
    bounds[bounds[,1] == -Inf, 1] <- 0
    bounds <- log(bounds)
    f <- function(x) fundamental_fitness(exp(trait_matrix(x, traits)), p)
  } else {
    f <- function(x) fundamental_fitness(trait_matrix(x, traits), p)
  }

  if (length(traits) == 1L) {
    if (!all(is.finite(bounds))) {
      stop("Starting value did not have finite fitness; finite bounds required")
    }
    ## The suppressWarnings here is for warnings like:
    ##
    ## Warning message:
    ## In optimise(f, interval = bounds, maximum = TRUE, tol = tol) :
    ##   NA/Inf replaced by maximum positive value
    ##
    ## which is probably the desired behaviour here.
    out <- suppressWarnings(optimise(f, interval=bounds, maximum=TRUE, tol=tol))
    ret <- out$maximum
    fitness <- out$objective
  } else {
    ## This is not very well tested, and the tolerance is not useful:
    out <- optim(rowMeans(bounds), f, method="L-BFGS-B",
                 lower=bounds[, "lower"], upper=bounds[, "upper"],
                 control=list(fnscale=-1, factr=1e10))
    ret <- out$par
    fitness <- out$value
  }

  if (log_scale) {
    ret <- exp(ret)
  }
  attr(ret, "fitness") <- fitness
  ret
}
