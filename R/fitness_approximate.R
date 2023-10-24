## The issue here is that we want to organise two different things:
## *points*, from which we can construct the approximate landscape and
## the actual approximate objects that can be looked up.  We'll store
## points within the object, but create actual functions only when
## needed.  The cost of doing this is cheap enough, even for the GPs.
## If that ends up being too costly, then I'll just memoise the
## generator function.

## The other problem is that we need to move a number of things back
## to calling functions so that subsequent calls are faster.

## The order of dependencies is:
##
##   node_schedule_times      } -- required for the fitness function
##   node_schedule_ode_times  } /
##   fitness_approximate_points   -- required for approximate 1d
##   fitness_approximate_slopes   -- required for approximate 2d
##
## So every time we call one of these functions we should push back
## any attributes of this type back.  This stuff is where references
## were definitely nice.
##
## Because the 1d and 2d functions return different bits of computed
## data it's particularly tricky.

## Same but for community:
community_prepare_fitness <- function(community) {
  if (!is.null(community$bounds)) {
    if (length(community) > 0L) {
      if (is.null(community$node_schedule_times)) {
        res <- build_schedule(community_parameters(community), ctrl = plant_control())
        community$node_schedule_times <- res$node_schedule_times
        community$node_schedule_ode_times <-
          res$node_schedule_ode_times
      } else if (is.null(community$node_schedule_ode_times)) {
        res <- run_scm(community_parameters(community))
        community$node_schedule_ode_times <- res$ode_times
      }
    }
  }
  community
}

community_prepare_approximate_fitness <- function(community) {
  community <- community_prepare_fitness(community)
  if (!is.null(community$bounds)) {
    if (length(community$trait_names) == 1L) { # 1d
      if (is.null(community$fitness_approximate_points)) {
        community$fitness_approximate_points <-
          fitness_approximate_points(community)
      }
    } else {
      if (length(community) > 0L) {
        community$fitness_approximate_slopes <-
          community_fitness_slopes(community)
      }
    }
  }
  community
}

community_check_fitness_prepared <- function(community, approximate) {
  n_spp <- length(community)
  n_traits <- length(community$trait_names)
  if (n_spp == 0) {
    return(TRUE)
  }
  ok_real <- (!is.null(community$node_schedule_times) &&
               !is.null(community$node_schedule_ode_times))
  ok_approximate <- (
    !approximate ||
    (n_traits == 1 && !is.null(community$fitness_approximate_points)) ||
    (n_traits > 1 &&  !is.null(community$fitness_2d_information)))

  ok_real && ok_approximate
}

community_assert_fitness_prepared <- function(community, approximate) {
  if (!community_check_fitness_prepared(community, approximate)) {
    stop("Unprepared")
  }
}

community_make_fitness <- function(community) {


  community_assert_fitness_prepared(community, FALSE)
  p <- community_parameters(community)
  if(is.null(community$hyperpar)) {
    hyperpar <- param_hyperpar(p)
  } else {
    hyperpar <- community$hyperpar
  }

  fitness <- function(x) {
    
    plant_log_assembler("Computing fitness landscape")

    plant::fitness_landscape(x, p, hyperpar = hyperpar, ctrl = plant_control())
  }
  fitness
}

## This one returns a function based on these points.  Eventually the
## stop branch will be replaced by a call to the generation function
## but I want to make sure we pass around the points as much as
## possible.
##' @export
community_fitness_approximate <- function(community) {
  community_assert_fitness_prepared(community, TRUE)

  pts <- community$fitness_approximate_points
  control <- fitness_approximate_control(community$fitness_approximate_control)
  if (is.null(community$fitness_approximate_points)) {
    stop("Points are not precomputed")
  }
  make_fitness_approximate_function(pts, control$type)
}

fitness_approximate_points <- function(community, bounds=NULL) {
  if (!inherits(community, "community")) {
    stop("Expected a community object")
  }
  if (is.null(bounds)) {
    bounds <- community$bounds
  }
  bounds <- check_bounds(bounds, finite=TRUE)
  if (nrow(bounds) != 1L) {
    stop("Only working for one trait at the moment")
  }
  
  control <- fitness_approximate_control(community$fitness_approximate_control)
  ## TODO: it's probably worth harvesting the ode times here?
  fitness <- community_make_fitness(community)
  if (control$type == "grid") {
    m <- fitness_approximate_points_grid(fitness, community, bounds, control)
  }
   else if (control$type == "gp") {
    m <- fitness_approximate_points_gp(fitness, community, bounds, control)
  } else {
    stop("Unknown approximate type ", control$type)
  }

  if (control$finite_only) {
    m <- m[is.finite(m[,2]),,drop=FALSE]
  }
  colnames(m) <- c(community$trait_names, "fitness")
  m
}

fitness_approximate_points_grid <- function(fitness, community, bounds, control) {

  if(length(community$trait_names) > 1) {
     stop("Fitness_approximate_points_grid only works for single trait:", community$trait_names)
  }

  x <- seq_log_range(bounds, control$n)

  if(!is.null(control$grid_include_residents) && control$grid_include_residents) {
    x <- sort(unique(c(x, community$traits)))
  }

  trait <- trait_matrix(x, community$trait_names)

  w <- fitness(trait)
  cbind(trait, w, deparse.level=0)
}

## TODO: why (here and above) is bounds not taken directly from the
## community?  Ditto control
fitness_approximate_points_gp <- function(fitness,
                                          community, bounds, control) {
  n_predict <- control$n_predict
  n_initial <- control$n_initial

  x_resident <- community$traits
  if (!is.null(control$x_seed)) {
    x_resident <- rbind(x_resident, matrix(control$x_seed, ncol=1L))
  }
  n_resident <- nrow(x_resident)


  objective <- function(x) {
    fitness(exp(x))
  }
  bounds <- log(bounds)
  x_resident <- log(x_resident)
  x <- trait_matrix(seq_range(bounds, n_predict),
                    community$trait_names)

  ## Set up the residents and any seed points as special data points
  ## to include in the search.  However, don't count this in the
  ## random set.  To avoid interface issues, these are computed
  ## *before* passing through to interp1d.
  src <- grail::data_source$new(rbind(x_resident, x), objective,
                                lower_limit=control$lower_limit,
                                verbose=FALSE)
  initial <- sample(n_predict, n_initial) + n_resident
  src$force_vals(c(seq_len(n_resident), initial))
  gp <- grail::interp1d(src,
                        n_initial + n_resident,
                        control$n + n_resident, control$n_each,
                        cost=control$cost)
  cbind(exp(t(gp$X)), t(gp$y), deparse.level=0)
}

##' @importFrom grader slope_info
community_fitness_slopes <- function(sys) {
  fitness <- community_make_fitness(sys)
  f_logspace <- function(x) {
    fitness(trait_matrix(exp(x), sys$trait_names))
  }
  apply(log(sys$traits), 1, function(x) grader::slope_info(f_logspace, x))
}

## Then, the next step is to generate the function to evaluate:
make_fitness_approximate_function <- function(pts, type) {
  trait <- pts[, 1, drop=FALSE]
  fitness <- pts[, 2, drop=FALSE]
  if (type == "grid") {
    splinefun_log(trait, fitness)
  } else {
    gp <- grail::gpreg("matern32")(log(trait), fitness)
    function(x) {
      gp$predict(log(cbind(x)))
    }
  }
}

##' List of parameters controlling
##' how approximate fitness landscapes are generated.
##'
##' @title Control how fitness landscapes are constructed
##' @param control List of values to change from default.
##' @return List with elements type, n, finite_only,
##' n_initial,, n_each, n_predict, lower_limit, cost,
##' x_seed.
##' @author Rich FitzJohn, Daniel Falster
##' @export
fitness_approximate_control <- function(control=NULL) {
  defaults <- list(type="grid",
                   n=50, # total number of points
                   finite_only=TRUE,
                   ## For gp:
                   n_initial=20,
                   n_each=5,
                   n_predict=500,
                   lower_limit=NULL,
                   cost=grail::cost_var_scaled_capped,
                   x_seed=NULL,
                   grid_include_residents=FALSE)
  modifyList(defaults, as.list(control))
}
