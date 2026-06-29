## Functions *related* to assembly, but not actually doing it.

#' Find evolutionary attractor for single species and trait.
#'
#' Find evolutionary attractor for single species and trait.
#' This is point at which selection gradient equals zero.
#' Currently solved using \code{uniroot}
#' @param community A \code{community} object to search within.
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
      community_add(trait_matrix(x, community$trait_names)) %>%
      community_demography() %>%
      community_selection_gradient()
    
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
#' using \code{community_demography}
#' @param community community object to use. 
#' @param dx Interval over which derivative is calculated
#' @param log_scale (currently disabled) Determines whether derivative is taken
#' with respect to raw or log-transformed x values. The latter
#' is useful when x is log-normally distributed.
#' @author Daniel Falster
#' @export
#' @return a community with selection gradient added.
community_selection_gradient <- function(community, dx=1e-04,
                                log_scale=TRUE) {

  msg <- sprintf("Calculating selection gradient for [%s] = [%s]",
    paste(community$trait_names, collapse = ", "), 
    paste(community$traits, collapse = ", ")
  )
  plant_log_assembler(msg)
  
  trait_names <- community$trait_names
  # get points needed for gradient
  points <- gradient_points(community$traits, d = dx, r = 1)
 
  # bind on current traits so we can return current fitness too
  xx <- 
    rbind(
      community$traits,
      points
    )
  
  # calculate fitness
  ff <- community$fitness_function(xx)

  # extract points for derivative
  y <- ff[-1]  
  dim(y) <- attr(points, "dim_y")
  # caluclate gradient using forward difference
  ret <- gradient_extrapolate(y, points)

  msg <- sprintf("Solved! Selection gradient for [%s] = [%s] is [%s]", 
    paste(community$trait_names, collapse = ", "), 
    paste(community$traits, collapse = ", "),
    paste(ret, collapse = ", ")
  )
  plant_log_assembler(msg)

  community[["resident_fitness"]] <- ff[1]
  community[["selection_gradient"]] <- ret
  community
}

## NOTE: max_fitness() and max_growth_rate() now live in
## R/community_fitness_solve_max.R, reimplemented on the community machinery.
## The previous plant-style versions here called the removed plant
## fitness_landscape_empty()/fundamental_fitness() and have been deleted.
