# Establish baseline parameters
plant_control <- function() {
  ctrl <- scm_base_control()
  #  ctrl$equilibrium_nsteps <- 20
  #  ctrl$equilibrium_solver_name <- "hybrid"
  ctrl$equilibrium_nsteps <- 20
  ctrl$equilibrium_solver_name <- "iteration"
  ctrl
}

f_assemble <- function(trait, disturbance_interval,
                       B_lf1,
                       filename,
                       community0 = NULL
                       ) {
  
  dir.create(dirname(filename), FALSE, TRUE)

  logfile <- gsub(".rds", ".txt", filename)
  plant_log_console(logfile)

  plant.assembly:::plant_log_assembler(
    paste0("starting assembly with disturbance_interval = ", disturbance_interval,
      "; B_lf1 = ", B_lf1,"; saving to ", filename))

  # Establish baseline parameters
  ctrl_plant <- scm_base_control()
  p <- scm_base_parameters("FF16")
  p$strategy_default$a_l1 <- 2.17
  p$strategy_default$a_l2 <- 0.5
  p$strategy_default$hmat <- 12
  p$strategy_default$lma <- 0.038

  
  p$max_patch_lifetime <- disturbance_interval * 3.5

  if(trait == "lma") {
    bounds0 <- bounds(lma = c(0.01, 1))
  } else {
    bounds0 <- bounds(hmat = c(0.5, 30))
  }

  # change some of the controls on assembly
  control <-
    assembler_control(
      list(
        run_type = "to_equilibrium",
        birth_type = "maximum",
        birth_move_tol = 1,
        compute_viable_fitness = FALSE,
        equilibrium_eps = ctrl_plant$equilibrium_eps
      )
    )

  # setup an empty community
  if (is.null(community0)) {
    community0 <-
      community_start(
        p,
        bounds0,
        hyperpar = make_FF16_hyperpar(B_lf1 = B_lf1),
        fitness_control = list(type = "gp")
      )
  }
  
  # setup the assembler
  obj <-
    assembler_start(
      community0,
      control = control,
      filename = filename
    )

  obj <- assembler_run(obj, nsteps = 20)

  saveRDS(obj, filename)

  TRUE
}

process_results <- function(pars) {
  results <-
    pars %>%
    mutate(
      label = paste(trait, disturbance_interval, B_lf1),
      obj = map(filename, readRDS),
      landscape = map(obj, ~ .x$community$fitness_points %>% as_tibble()),
      res = map2(obj, landscape, ~ tibble(trait_value = .x$community$traits[, 1], fitness = approx(.y[[1]], .y[[2]], .x$community$traits[, 1])$y))
    )

  landscapes <-
    results %>%
    select(label, trait, disturbance_interval, B_lf1, landscape) %>%
    unnest(landscape) %>%
    mutate(trait_value = ifelse(is.na(lma), hmat, lma))

  residents <-
    results %>%
    select(label, trait, disturbance_interval, B_lf1, res) %>%
    unnest(res)

  list(
    results = results,
    landscapes = landscapes,
    residents = residents
  )
}


plot_landscapes <- function(landscapes, residents) {
  ggplot(landscapes, aes(trait_value, fitness)) +
    geom_line() +
    theme_classic() +    
    scale_x_log10() +
    geom_point(data = residents, col = "red")
}

plot_community_landscape <- function(community) {
  
  landscapes <- community$fitness_points %>% as_tibble() %>%
    mutate(trait_value = ifelse(is.na(lma), hmat, lma))

  residents <- tibble(trait_value = community$traits[, 1], fitness = approx(landscapes[[1]], landscapes[[2]], community$traits[, 1])$y)
  
  plot_landscapes(landscapes, residents)
}



plot_landscape_i <- 
  function(i, landscapes, residents) {
  ggplot(landscapes |> filter(step==i), 
      aes(trait_value, fitness)) +
    geom_line() +
    theme_classic() +    
    scale_x_log10(limits=c(0.009,1.02)) +
    geom_hline(yintercept = 0, lty="dashed") +
    geom_point(data = residents %>% filter(step==i), col="red") + theme_classic() + 
  ylim(c(-50,15)) +
  xlab("LMA (kg/m2)") +
  ylab("Fitness")
}
