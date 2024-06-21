plant_log_viable <- function(...) {
  plant_log_info(..., routine = "viable")
}

plant_log_inviable <- function(...) {
  plant_log_info(..., routine = "inviable")
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
community_viable_fitness_1D <- function(bounds, p, x=NULL, log_scale=TRUE, dx=1) {
  bounds <- check_bounds(bounds)
  traits <- rownames(bounds)
  
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
    stop("review implementation of 2D viable fitness landscape")
    #out <- positive_2d(f, x, lower=bounds[,1], upper=bounds[,2])
    #ret <- t(apply(out$points[out$labels,], 2, range))
  }
  if (log_scale) {
    ret <- exp(ret)
  }
  colnames(ret) <- c("lower", "upper")
  rownames(ret) <- traits
  ret
}
