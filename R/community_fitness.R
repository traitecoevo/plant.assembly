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
fitness_landscape <- function(trait_matrix, p, hyperpar=param_hyperpar(p), log_fitness=TRUE, ctrl = scm_base_control()) {
  
  n_residents <- length(p$strategies)
  mutant_birth_rates <- rep(0, nrow(trait_matrix))
  p_mutants <- mutant_parameters(trait_matrix, p, hyperpar,
                                 birth_rate_list = mutant_birth_rates)

  if (n_residents > 0L) {
    # if there's a resident, use the saved environment to calculate mutant fitness
    ctrl$save_RK45_cache = T  

    scm <- run_scm(p, use_ode_times=length(p$node_schedule_ode_times) > 0, ctrl = ctrl)
    scm$run_mutant(p_mutants)
    net_reproduction_ratios <- scm$net_reproduction_ratios
  } else {
    # otherwise just run mutants with zero birth rate
    scm <- run_scm(p_mutants, use_ode_times=length(p$node_schedule_ode_times) > 0,
                   ctrl = ctrl)
    net_reproduction_ratios <- scm$net_reproduction_ratios
  }

  if (log_fitness) {
    log(net_reproduction_ratios)
  } else {
    net_reproduction_ratios
  }
}
