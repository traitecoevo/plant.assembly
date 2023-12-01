
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
##' @param trait_matrix A matrix of traits corresponding to mutants to
##' introduce into the light environment constructed by the residents
##' in \code{p}.
##' @param p Parameters object.  Needs to contain residents with their
##' incoming seed rain.
##' @param hyperpar Hyperparameter function to use. By default links
##' to standard function for this strategy type.
##' @param log_fitness Logical; if \code{TRUE} report per capita
##' seed rain rather than fitness.
##' @param ctrl A plant control object
##' @return Vector with the output seed rain.  Mutants have an
##' arbitrary seed rain of one, so this is the rate of seed
##' production per capita.
##' @author Rich FitzJohn
##' @export
plant_community_fitness_landscape <- function(community, bounds = community$bounds, n = community$fitness_control$n) {

  plant_log_assembler(sprintf(
    "Calulcating fitness landscape for %d strategy communtiy, %d points", nrow(community$traits), n))

  x <- seq_log_range(bounds, n)

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
    # otherwise just run mutants with zero birth rate
    scm <- run_scm(p_mutants, ctrl = ctrl,
      use_ode_times = length(p$node_schedule_ode_times) > 0)
    community$resident_fitness <- numeric()
  }

  community$fitness_function <- 
    function(x) {
      traits <- trait_matrix(x, community$trait_names)

      p_mutants <- mutant_parameters(traits, p, hyperpar,
        birth_rate_list = rep(0, nrow(traits)) )

      scm$run_mutant(p_mutants)
      log(scm$net_reproduction_ratios)
  }

  community
}



# connect functions
community_parameters <- plant_community_parameters
community_run <- plant_community_run
community_run_to_equilibrium <- plant_community_run_to_equilibrium

community_fitness_landscape <- plant_community_fitness_landscape
community_viable_bounds <- plant_community_viable_bounds
