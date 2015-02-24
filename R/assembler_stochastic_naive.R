make_births_stochastic_naive <- function(control) {
  vcv <- control$vcv
  if (is.null(vcv) || !is.matrix(vcv)) {
    stop("A vcv must be provided")
  }
  mutation    <- make_mutation_stochastic_naive(control$n_mutants, vcv)
  immigration <- make_immigration_stochastic_naive(control$n_immigrants)
  function(sys) {
    to_add <- rbind(mutation(sys), immigration(sys))
    if (control$check_positive && nrow(to_add) > 0) {
      fitness <- community_make_fitness(test)
      w <- fitness(to_add)
      keep <- w >= 0.0
      to_add <- to_add[keep, , drop=FALSE]
      attr(to_add, "fitness") <- w[keep]
    }
    to_add
  }
}

make_deaths_stochastic_naive <- function(control) {
  eps <- control$dead_seed_rain
  check_inviable <- control$check_inviable
  function(sys) {
    sys <- community_drop(sys, sys$seed_rain < eps)
    if (check_inviable) {
      sys <- community_drop_inviable(sys)
    }
    sys
  }
}

## Mutation: Draw (on average) n_mutants from the population with a
## mutational variance of vcv on the log scale.
make_mutation_stochastic_naive <- function(n_mutants, vcv) {
  n_traits <- ncol(vcv)
  blank <- matrix(nrow=0, ncol=n_traits)
  function(sys) {
    if (length(sys) == 0) {
      return(blank)
    }
    n <- rpois(1, n_mutants)
    if (n == 0) {
      return(blank)
    }
    traits <- sys$traits
    weights <- sys$seed_rain
    i <- sample(length(weights), n, replace=TRUE, prob=weights)
    unname(exp(log(traits[i,,drop=FALSE]) + mvtnorm::rmvnorm(n, sigma=vcv)))
  }
}

## Immigration: Same from the bounds, in log space.
make_immigration_stochastic_naive <- function(n_immigrants) {
  function(sys) {
    bounds <- sys$bounds
    n_traits <- length(sys$trait_names)
    if (is.null(bounds)) {
      ret <- matrix(nrow=0, ncol=n_traits)
    } else {
      lower <- log(bounds[,1])
      range <- log(bounds[,2]) - lower
      n <- rpois(1, n_immigrants)
      if (n == 0) {
        ret <- matrix(nrow=0, ncol=n_traits)
      } else {
        u <- t(lhs::randomLHS(n, n_traits))
        ret <- exp(t(lower + range * u))
      }
    }
    colnames(ret) <- sys$trait_names
    ret
  }
}

##' Generate a mutational variance covariance matrix where each trait
##' mutates by a fraction of its total range, and where there is no
##' covariance among traits.  This is a helper function for use with
##' \code{\link{assembler_stochastic_naive}}
##' @title Generate mutational covariance matrices
##' @param x A community object, or appropriate matrix of bounds
##' @param p Fraction of trait range to mutate (vectors will be
##' recycled according R's usual rules).
##' @export
mutational_vcv_proportion <- function(x, p=0.001) {
  if (inherits(x, "community")) {
    x <- x$bounds
  }
  bounds <- check_bounds(x)
  if (!all(is.finite(bounds))) {
    stop("All bounds must be finite")
  }
  p * diag(nrow(bounds)) * as.numeric(diff(t(log(bounds))))
}
