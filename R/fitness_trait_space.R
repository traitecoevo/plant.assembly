## These functions involve "trait space" as a concept, and as such
## care about hyperparameters; ways of converting things like lma into
## a set of parameters.  This is handled by strategy_list, which all
## of these functions use.

##' Find point of maximum fitness in empty fitness landscape within a
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

##' Compute region of positive fitness.  This will be the values where
##' fitness is approximately zero.
##'
##' @title Compute Region of Positive Fitness
##' @param bounds Matrix of bounds; two columns corresponding to the
##' lower and upper limits, each row corresponds to a trait (the name
##' will be used).
##' @param p Parameters object to use.  Importantly, the
##' \code{strategy_default} element gets used here.
##' @param x Initial value - If not given, then the value from the
##' default strategy within \code{p} is used.
##' @param log_scale Is the parameter naturally on a log scale?  If
##' so, this will greatly speed things up.  Can be a vector of length
##' `\code{nrow(bounds)}`
##' @param dx Amount to step the trait.  If \code{log_scale} is
##' \code{TRUE}, this is on a log scale.
##' @export
##' @author Rich FitzJohn
viable_fitness <- function(bounds, p, x=NULL, log_scale=TRUE, dx=1) {
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)
  if (is.null(x)) {
    x <- unlist(p$strategy_default[traits])
  }
  x <- check_point(x, bounds)
  w <- fundamental_fitness(x, p)

  if (w < 0) {
    plant_log_viable("Starting value had negative fitness, looking for max")
    x <- max_fitness(bounds, p, log_scale)
    w <- attr(x, "fitness")
    plant_log_viable(sprintf("\t...found max fitness at %s (w=%2.5f)",
                             paste(formatC(x), collapse=", "), w))
    if (w < 0) {
      return(NULL)
    }
  }

  if (log_scale) {
    bounds[bounds[,1] == -Inf, 1] <- 0
    bounds <- log(bounds)
    x <- log(x)
    f <- function(x) {
      fundamental_fitness(exp(trait_matrix(x, traits)), p)
    }
  } else {
    f <- function(x) {
      fundamental_fitness(trait_matrix(x, traits), p)
    }
  }

  if (length(traits) == 1L) {
    out <- positive_1d(f, x, dx, lower=bounds[,1], upper=bounds[,2])
    ret <- rbind(out, deparse.level=0)
  } else {
    ## TODO: It seems a crying shame to throws all this information
    ## away.  I'd like to create some sort of polygon around the
    ## positive points but need to get the bilinear interpolation
    ## done.  Until then, this is lost.
    
    stop("review implementation of 2D landscape")
    ##out <- positive_2d(f, x, lower=bounds[,1], upper=bounds[,2])
    ret <- t(apply(out$points[out$labels,], 2, range))
  }
  if (log_scale) {
    ret <- exp(ret)
  }
  colnames(ret) <- c("lower", "upper")
  rownames(ret) <- traits
  ret
}


##' @importFrom stats uniroot
positive_1d <- function(f, x, dx, lower = -Inf, upper = Inf, tol = 1e-3) {
  root <- function(b, type) {
    x <- b[[type]]$x
    fx <- b[[type]]$fx
    if (prod(fx[1:2]) < 0) {
      ## The suppressWarnings is here to stop the warning
      ##   -Inf replaced by maximally negative value
      ## which we're actually OK with.
      suppressWarnings(uniroot(f, x,
        f.lower = fx[[1]], f.upper = fx[[2]],
        tol = tol
      )$root)
    } else {
      if (type == "lower") x[[1]] else x[[2]]
    }
  }
  b <- positive_1d_bracket(f, x, dx, lower, upper)
  c(root(b, "lower"), root(b, "upper"))
}

positive_1d_bracket <- function(f, x, dx, lower, upper, grow = 2) {
  fx <- f(x)
  if (fx < 0) {
    stop("Don't yet support doing this with no positive values")
  }

  L <- U <- x
  fL <- fU <- fx
  dL <- dU <- dx

  bracket <- function(x, dx, bound) {
    cleanup <- function(x, x_next, fx, fx_next) {
      if (dx < 0) {
        x <- c(x_next, x)
        fx <- c(fx_next, fx)
      } else {
        x <- c(x, x_next)
        fx <- c(fx, fx_next)
      }
      list(x = x, fx = fx)
    }
    hit_bounds <- FALSE
    repeat {
      x_next <- x + dx
      if ((dx < 0 && x_next < bound) || (dx > 0 && x_next > bound)) {
        x_next <- bound
        hit_bounds <- TRUE
      }
      fx_next <- f(x_next)
      if (fx_next < 0 || hit_bounds) {
        return(cleanup(x, x_next, fx, fx_next))
      } else {
        x <- x_next
        fx <- fx_next
        dx <- dx * grow
      }
    }
  }

  list(
    lower = bracket(x, -dx, lower),
    upper = bracket(x, dx, upper)
  )
}

