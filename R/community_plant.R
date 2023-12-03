
plant_default_assembly_control <- function(equilibrium_nsteps = 100,equilibrium_solver_name = "iteration") {
  plant_control <- plant::scm_base_control()  
  plant_control$equilibrium_nsteps <- equilibrium_nsteps
  plant_control$equilibrium_solver_name <- "iteration"
  plant_control
}

plant_default_assembly_pars <- function(hmat = 10, max_patch_lifetime = 60, fixed_RA = FALSE) {
  p <- scm_base_parameters("FF16")
  p$strategy_default$a_l1 <- 2.17
  p$strategy_default$a_l2 <- 0.5
  p$strategy_default$hmat <- hmat
  p$max_patch_lifetime <- max_patch_lifetime
  if(fixed_RA) {
    # Fixed allocation to reproduction
    p$strategy_default$a_f1 = 0.5
    p$strategy_default$a_f2 = 0
  }
  p
}

plant_community_hyperpar <- function(community) {
 if (is.null(community$model_support$hyperpar)) {
   hyperpar <- plant::param_hyperpar(community$model_support$p)
 } else {
   hyperpar <- community$model_support$hyperpar
 }
 hyperpar
}


plant_community_parameters <- function(community) {
  
  p <- community$model_support$p

  hyperpar <- plant_community_hyperpar(community)

  if(nrow(community$traits) > 0)
    p$strategies <- plant::strategy_list(community$traits, p, birth_rate_list = community$birth_rate, hyperpar = hyperpar)

  if (!is.null(community$model_support$node_schedule_times)) {
    p$node_schedule_times <- community$model_support$node_schedule_times
  }
  if (!is.null(community$model_support$node_schedule_ode_times)) {
    p$node_schedule_ode_times <- community$model_support$node_schedule_ode_times
  }

  p
}

plant_community_run <- function(community) {
  if (length(community) > 0L) {

    p <- plant:::build_schedule(plant_community_parameters(community), ctrl = community$model_support$plant_control)
    community$birth_rate <- attr(p, "offspring_production")

    community$model_support$node_schedule_times <- p$node_schedule_times
    community$model_support$node_schedule_ode_times <- p$node_schedule_ode_times
    community$fitness_points <- NULL
  }

  plant_community_update_fitness_function(community)
}

plant_community_run_to_equilibrium <- function(community) {
  if (length(community) > 0L) {
    p <- equilibrium_birth_rate(plant_community_parameters(community), 
            ctrl = community$model_support$plant_control)

    community$birth_rate <- attr(p, "offspring_production")
    community$model_support$node_schedule_times <- p$node_schedule_times
    community$model_support$node_schedule_ode_times <- p$node_schedule_ode_times
    community$fitness_points <- NULL
  }

  plant_community_update_fitness_function(community)
}

plant_community_viable_bounds <- function(community) {
  if (length(community) > 0) {
    stop("You don't want to run this on an existing community")
  }
  community$bounds <- viable_fitness(community$bounds, plant_community_parameters(community))
  community
}


##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param community A community object
##' @param bounds a bounds object
##' @param npts number of points 
##' @author Daniel Falster
##' @export
community_fitness_landscape <- function(community, bounds = community$bounds, npts = community$fitness_control$n) {

  if(is.null(community$fitness_function)) {
    community <- community %>% community_run()
  }
  plant_log_assembler(sprintf(
    "Calulcating fitness landscape for %d strategy communtiy, %d points", nrow(community$traits), npts))

  x <- seq_log_range(bounds, npts)

  # add residents
  x <- sort(unique(c(x, community$traits)))

  y <- community$fitness_function(x)

  community$fitness_points <- 
    dplyr::tibble(x = x, fitness = y) %>% 
    mutate(resident = ifelse(x %in% community$traits, TRUE, FALSE))
  
  names(community$fitness_points)[1] <- community$trait_names
 
  community
}

plant_community_update_fitness_function <- function(community) {
  
  p <- plant_community_parameters(community)
  hyperpar <- plant_community_hyperpar(community)

  ctrl <- community$model_support$plant_control
  ctrl$save_RK45_cache <- T
  
  if (length(p$strategies) > 0L) {
    # if there's a resident, use the saved environment to calculate mutant fitness
    scm <- run_scm(p, ctrl = ctrl,
      use_ode_times = length(p$node_schedule_ode_times) > 0)
    community$resident_fitness <- log(scm$net_reproduction_ratios)
    } else {
    community$resident_fitness <- numeric()
  }

  community$fitness_function <- 
    function(x) {
      traits <- trait_matrix(x, community$trait_names)

      p_mutants <- mutant_parameters(traits, p, hyperpar,
        birth_rate_list = rep(0, nrow(traits)))
      
      if (length(p$strategies) > 0L) {
        scm$run_mutant(p_mutants)
        ret <- log(scm$net_reproduction_ratios)
      } else {
        # otherwise just run mutants with zero birth rate
        scm <- run_scm(p_mutants, ctrl = ctrl, use_ode_times = 0)
        ret <- log(scm$net_reproduction_ratios)
      }
      
    ret
  }

  community
}

##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param community A community object
##' @export
plant_community_check_for_inviable_strategies <- function(community) {
  # TODO: Remove direct depoency o n plant, rewrite to use existing community fitness functions

  p <- community_parameters(community)
  ctrl <- community$model_support$plant_control

  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$birth_rate)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- ctrl$equilibrium_extinct_birth_rate
  ## TODO: This was ctrl$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  birth_rate <- sapply(p$strategies, function(s) s$birth_rate_y, simplify = TRUE)
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- plant:::make_equilibrium_runner(p, ctrl = ctrl)
  offspring_production <- runner(birth_rate)

  test <- which(offspring_production < birth_rate &
    birth_rate < max(offspring_production) * eps_test)
  test <- test[order(offspring_production[test])]

  drop <- logical(length(offspring_production))

  for (i in test) {
    plant:::plant_log_inviable(paste("Testing species", i),
      stage = "testing", species = i
    )
    x <- offspring_production
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant:::plant_log_inviable(paste("Removing species", i),
        stage = "removing", species = i
      )
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      offspring_production <- res
    }
  }

  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | offspring_production < eps

  attr(offspring_production, "drop") <- drop
  offspring_production
}

# connect functions
community_parameters <- plant_community_parameters
community_run <- plant_community_run
community_run_to_equilibrium <- plant_community_run_to_equilibrium
community_viable_bounds <- plant_community_viable_bounds
community_check_for_inviable_strategies <- plant_community_check_for_inviable_strategies

