community_new_types_stochastic <- function(sys, control) {
  plant_log_assembler("Adding new types")

  vcv <- control$vcv
  if (is.null(vcv) || !is.matrix(vcv)) {
    stop("A vcv must be provided")
  }
  mutation    <- make_mutation_stochastic_naive(control$n_mutants, vcv)
  immigration <- make_immigration_stochastic_naive(control$n_immigrants)

  
  to_add <- rbind(mutation(sys), immigration(sys))
  if (control$check_positive && nrow(to_add) > 0) {
    w <- sys$fitness_function(to_add)#fitness(to_add)
    keep <- w >= 0.0
    to_add <- to_add[keep, , drop=FALSE]
    attr(to_add, "fitness") <- w[keep]
  }
  to_add
}

## Mutation: Draw (on average) n_mutants from the population with a
## mutational variance of vcv, applied on the community's trait scale
## (log for strictly-positive traits, linear when trait_scale = "linear";
## see community_trait_transform()).
make_mutation_stochastic_naive <- function(n_mutants, vcv) {
  n_traits <- ncol(vcv)
  blank <- matrix(nrow=0, ncol=n_traits)
  function(sys) {
    if (length(sys) == 0) {
      return(blank)
    }
    n_mutants_actual <- rpois(1, n_mutants)
    if (n_mutants_actual == 0) {
      return(blank)
    }
    tf <- community_trait_transform(sys)
    traits <- sys$traits
    weights <- sys$birth_rate
    i <- sample(length(weights), n_mutants_actual, replace = TRUE, prob = weights)
    unname(tf$inv(tf$fwd(traits[i, , drop = FALSE]) +
                    mvtnorm::rmvnorm(n_mutants_actual, sigma = vcv)))
  }
}

## Immigration: draw uniformly across the bounds, on the community's trait
## scale (log or linear; see community_trait_transform()).
make_immigration_stochastic_naive <- function(n_immigrants) {
  function(sys) {
    bounds <- sys$bounds
    n_traits <- length(sys$trait_names)
    if (is.null(bounds)) {
      ret <- matrix(nrow=0, ncol=n_traits)
    } else {
      tf <- community_trait_transform(sys)
      lower <- tf$fwd(bounds[,1])
      range <- tf$fwd(bounds[,2]) - lower
      n <- rpois(1, n_immigrants)
      if (n == 0) {
        ret <- matrix(nrow=0, ncol=n_traits)
      } else {
        u <- t(lhs::randomLHS(n, n_traits))
        ret <- tf$inv(t(lower + range * u))
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
##' @details The trait range is measured on the community's trait scale: the
##' log scale by default (so the variance is a proportion of the log range, as
##' for strictly-positive traits), or the linear scale when the community was
##' created with \code{trait_scale = "linear"}. A bare bounds matrix has no
##' scale attached and is treated as log, preserving the historical behaviour.
##' @export
mutational_vcv_proportion <- function(x, p=0.001) {
  fwd <- log                                   # default scale for bare bounds
  if (inherits(x, "community")) {
    fwd <- community_trait_transform(x)$fwd
    x <- x$bounds
  }
  bounds <- check_bounds(x)
  if (!all(is.finite(bounds))) {
    stop("All bounds must be finite")
  }
  p * diag(nrow(bounds)) * as.numeric(diff(t(fwd(bounds))))
}
