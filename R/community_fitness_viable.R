plant_log_viable <- function(...) {
  plant_log_info(..., routine = "viable")
}

plant_log_inviable <- function(...) {
  plant_log_info(..., routine = "inviable")
}

##' Compute the region of positive (fundamental) fitness for a community.
##'
##' Finds the trait range over which a strategy has positive invasion fitness
##' into the community's environment. This is normally run on an *empty*
##' community, in which case it is the fundamental niche: the range where a lone
##' strategy can persist with no competitors.
##'
##' Reimplemented on the community machinery: it uses
##' \code{community$fitness_function} (built by
##' \code{plant_community_update_fitness_function}) as the fitness function,
##' rather than the removed plant \code{fundamental_fitness()}/\code{viable_fitness()}.
##'
##' @title Compute region of positive fitness for a community
##' @param community A \code{community} object (usually empty).
##' @param x Initial trait value. If \code{NULL}, the value from the default
##'   strategy in \code{community$model_support$p} is used.
##' @param log_scale Is the trait naturally on a log scale? If so this speeds
##'   up the search considerably.
##' @param dx Amount to step the trait when bracketing (on the log scale when
##'   \code{log_scale} is \code{TRUE}).
##' @return A bounds matrix (\code{lower}/\code{upper} columns, one row per
##'   trait), or \code{NULL} if no positive-fitness region was found.
##' @export
##' @author Rich FitzJohn, Daniel Falster
community_viable_fitness_1D <- function(community, x = NULL,
                                        log_scale = TRUE, dx = 1) {
  bounds <- check_bounds(community$bounds)
  traits <- rownames(bounds)
  if (length(traits) != 1L) {
    stop("review implementation of 2D viable fitness landscape")
  }

  ## Fundamental fitness == invasion fitness into this (empty) community.
  if (is.null(community$fitness_function)) {
    community <- plant_community_update_fitness_function(community)
  }
  fitness <- function(trait_value) {
    community$fitness_function(as.numeric(trait_value))
  }

  ## Default starting point: the default strategy's trait value.
  if (is.null(x)) {
    x <- unlist(community$model_support$p$strategy_default[traits])
  }
  x <- check_point(x, bounds)
  w <- fitness(x)

  if (w < 0) {
    plant_log_viable("Starting value had negative fitness, looking for max")
    x <- community_max_fitness_1D(fitness, bounds, log_scale)
    w <- attr(x, "fitness")
    plant_log_viable(sprintf("\t...found max fitness at %s (w=%2.5f)",
                             paste(formatC(x), collapse = ", "), w))
    if (w < 0) {
      return(NULL)
    }
  }

  if (log_scale) {
    bounds[bounds[, 1] == -Inf, 1] <- 0
    bounds <- log(bounds)
    x <- log(as.numeric(x))
    f <- function(z) fitness(exp(z))
  } else {
    x <- as.numeric(x)
    f <- function(z) fitness(z)
  }

  out <- positive_1d(f, x, dx, lower = bounds[, 1], upper = bounds[, 2])
  ret <- rbind(out, deparse.level = 0)

  if (log_scale) {
    ret <- exp(ret)
  }
  colnames(ret) <- c("lower", "upper")
  rownames(ret) <- traits
  ret
}

## 1D maximisation of a fitness function over `bounds`, used when the default
## starting point has negative fitness. Operates directly on the community
## fitness function (no dependency on the removed plant max_fitness()).
community_max_fitness_1D <- function(fitness, bounds, log_scale = TRUE) {
  if (log_scale) {
    bounds[bounds[, 1] == -Inf, 1] <- 0
    b <- log(bounds)
    f <- function(z) fitness(exp(z))
  } else {
    b <- bounds
    f <- function(z) fitness(z)
  }
  if (!all(is.finite(b))) {
    stop("Starting value did not have finite fitness; finite bounds required")
  }
  ## suppressWarnings: optimise warns "NA/Inf replaced by maximum positive
  ## value" for inviable trait values, which is the behaviour we want.
  out <- suppressWarnings(optimise(f, interval = b, maximum = TRUE, tol = 1e-3))
  ret <- if (log_scale) exp(out$maximum) else out$maximum
  attr(ret, "fitness") <- out$objective
  ret
}

## --- pure numerical root helpers (reimplemented, no plant dependency) -------
## Find the interval around `x` (which must have f(x) >= 0) where f crosses zero,
## i.e. the bounds of the positive-fitness region. `positive_1d_bracket` brackets
## the sign change in each direction; `positive_1d` then refines with uniroot.

positive_1d <- function(f, x, dx, lower = -Inf, upper = Inf, tol = 1e-3) {
  root <- function(b, type) {
    xs <- b[[type]]$x
    fx <- b[[type]]$fx
    if (prod(fx[1:2]) < 0) {
      ## suppressWarnings: -Inf replaced by maximally negative value, which is OK.
      suppressWarnings(uniroot(f, xs,
                               f.lower = fx[[1]], f.upper = fx[[2]],
                               tol = tol)$root)
    } else {
      if (type == "lower") xs[[1]] else xs[[2]]
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

  list(lower = bracket(x, -dx, lower),
       upper = bracket(x, dx, upper))
}
