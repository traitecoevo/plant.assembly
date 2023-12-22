##' Construct a fitness landscape.
##'
##' @title Fitness Landscape
##' @param community A community object
##' @param method used to construct landscape
##' @param bounds a bounds object
##' @param npts number of points
##' @author Daniel Falster
##' @rdname community_fitness_landscape
##' @export
community_fitness_landscape <- function(community, method = community$fitness_control$method, ...) {

  if (is.null(community$fitness_function)) {
    community <- community %>% community_run()
  }
  plant_log_assembler(sprintf(
    "Calulcating fitness landscape for %d strategy communtiy using %s", nrow(community$traits), method))  
  
  if(method == "grid") {
    community_fitness_landscape_grid(community, ...)
  } else if(method == "bayes") {
    community_fitness_landscape_bayes(community, ...)
  } else {
    stop("Unknown fitness landscape method")
  }
}

community_fitness_landscape_grid <- function(community, bounds = community$bounds, npts = community$fitness_control$npts) {

  x <- seq_log_range(bounds, npts)

  # add residents
  x <- sort(unique(c(x, community$traits)))

  y <- community$fitness_function(x)

  community$fitness_points <-
    dplyr::tibble(x = x, fitness = y) %>%
    dplyr::mutate(resident = ifelse(x %in% community$traits, TRUE, FALSE))

  names(community$fitness_points)[1] <- community$trait_names

  community
}

##' Construct a fitness landscape using Bayesian optimisation.
## Uses the mlr3mbo package

community_fitness_landscape_bayes <- function(community, bounds = community$bounds, npts = community$fitness_control$npts, ninit = community$fitness_control$ninit) {
 
  set.seed(1)
  
  require(mlr3learners) # need this to make "regr.km" learner accessible

  obfun <- bbotk::ObjectiveRFun$new(
    fun = function(xs) list(community$fitness_function(exp(xs$x))),
    domain = paradox::ps(x = paradox::p_dbl(lower = log(community$bounds[1]), upper = log(community$bounds[2]))),
    codomain = paradox::ps(y = paradox::p_dbl(tags = "maximize"))
  )

  surrogate <- mlr3mbo::srlrn(mlr3::lrn("regr.km", covtype = "matern3_2", control = list(trace = FALSE)))

  optimizer <- bbotk::opt("mbo",
    loop_function = mlr3mbo::bayesopt_ego,
    surrogate = surrogate,
    acq_function = mlr3mbo::acqf("ei"),
    acq_optimizer = mlr3mbo::acqo(
      bbotk::opt("nloptr", algorithm = "NLOPT_GN_ORIG_DIRECT"),
      terminator = bbotk::trm("stagnation", iters = 100, threshold = 1e-5)
    )
  )

  instance <- bbotk::OptimInstanceSingleCrit$new(
    objective = obfun,
    terminator = bbotk::trm("evals", n_evals = npts)
  )

  # Initial data -- 
  # space npts and add residents
  x <- sort(unique(c(seq_log_range(bounds, min(ninit, npts)), community$traits)))

  initial_design <- data.table::data.table(x = log(x))
  instance$eval_batch(initial_design)

  ##input existing knowledge
  #initial_design <- data.table(x = log(community$fitness_points[[1]]), y = community$fitness_points[["fitness"]], batch_nr = 1)
  # instance$archive$data <- initial_design

  # run optimisation
  optimizer$optimize(instance)

  # Store points
  community$bayes_archive <- instance$archive

  community$fitness_points <-
    instance$archive$data %>% 
    dplyr::as_tibble() %>%
    dplyr::mutate(x = exp(x)) %>%  #back transform x
    dplyr::select(x = x, fitness = y, batch_nr) %>%
    dplyr::arrange(x) %>%
    dplyr::mutate(resident = ifelse(x %in% community$traits, TRUE, FALSE))
  names(community$fitness_points)[1] <- community$trait_names

  # Store surrogate
  community <- community_fitness_create_surrogate(community)

  #community$fitness_surrogate(0.01)
  #community$fitness_function(0.01)
  
  # # make predictions
  # xdt <- data.table::data.table(x = seq_range(log(community$bounds), length.out = 101))
  # surrogate_pred <- surrogate$predict(xdt)
  # pred <- bind_cols(xdt %>% as_tibble(), surrogate_pred %>% as_tibble()) %>%
  #   rename(y = mean)

  # # plot true function in black
  # # surrogate prediction (mean +- se in grey)
  # # known optimum in darkred
  # # found optimum in darkgreen
  # ggplot(aes(x, y), data = pred) +
  #   geom_point(data = instance$archive$data %>% as_tibble()) +
  #   geom_line(col = "red") +
  #   geom_ribbon(aes(x = x, ymin = y - se, ymax = y + se), fill = "grey", alpha = 0.2) +
  #   theme_minimal()
  community
}

community_fitness_create_surrogate <- function(community) {

  community$fitness_surrogate_obj <- 
    mlr3mbo::srlrn(mlr3::lrn("regr.km", covtype = "matern3_2", control = list(trace = FALSE)), archive = community$bayes_archive)

  community$fitness_surrogate_obj$update()

    
  community$fitness_surrogate <- function(x, se = FALSE) {
    xdt <- data.table::data.table(x = log(x))

    surrogate_pred <-
      community$fitness_surrogate_obj$predict(xdt) %>%
      dplyr::as_tibble() %>%
      dplyr::rename(y = mean)

    if (!se) {
      surrogate_pred <- surrogate_pred[["y"]]
    }

    surrogate_pred
  }
  
  community
}
