## Define models to run and parameters to vary
model_info <- function(name, res=80, default_pars = list(a_l1 = 2.17, a_l2 = 0.5)) {

  dat <- list()

  dat$narea <- list(trait_names="narea",
                  names=c("disturbance_mean_interval", "B_lf5"),
                  pars=default_pars)

  dat$lma <- list(trait_names="lma",
                  names=c("disturbance_mean_interval", "B_kl2"),
                  pars=default_pars)

  dat$rho <- list(trait_names="rho",
                  names=c("disturbance_mean_interval", "B_ks2"),
                  pars=default_pars)

  dat$rho_s <- list(trait_names="rho",
                  names=c("disturbance_mean_interval", "B_dI2"),
                  pars=default_pars)

  dat$omega <- list(trait_names="omega",
                   names=c("disturbance_mean_interval", "B_f2"),
                   pars=default_pars)

  if (name %in% names(dat)) {
    ret <- modifyList(dat[[name]], list(model=name))
  } else {
    stop(sprintf("model %s not found", name))
  }

  ## If we want more than disturbance_mean_interval to be logged, then
  ## change this to allow that.
  ret$log <- ret$names == "disturbance_mean_interval"
  ret$bounds <- trait_bounds(ret$trait_names)
  ## Fixed for 1d:
  ret$x <- plant::seq_log_range(ret$bounds, res)
  ret
}


## This is now nicely generic; we just need a little bit of data in
## `info` and we can construct everything required for run_model.
run_model_info <- function(pars, info,...) {
  if (length(pars) != length(info$names)) {
    stop(sprintf("Expected %d parameters", length(info$names)))
  }
  run_model(info$trait_names,
            c(setNames(as.list(pars), info$names), info$pars),
            ...)
}

run_model <- function(trait_names, pars, path, nsteps = 20L) {

  ## 1: logging
  filenames <- list(
    log = paste0(path, ".log"),
    rds = paste0(path, ".rds"),
    pars = paste0(path, "_pars.rds")
    )

  if(file.exists(filenames$log)){
    message(sprintf("gah, %s already exists, aborting", filenames$log))
    return()
  }

  
  dir.create(dirname(filenames$log), FALSE, TRUE)
  saveRDS(pars, filenames$pars)
  loggr::log_file(filenames$log)
  plant::plant_log_console()
  on.exit(stop_logging())

  p <- assembly_parameters(pars=pars, type = "FF16")
  
  p$hyperpar = make_new_hyperpar()
  log <- plant:::plant_log_info
  log(sprintf("Saving results to %s", path))
  
  log(sprintf("Assembler for traits %s", paste(trait_names, collapse=", ")))
  log(sprintf("\tmodifying parameters:\n%s",
              paste(sprintf("\t- %s: %s", names(pars), unlist(pars)),
                    collapse="\n")))
  log(sprintf("\trunning %s steps", nsteps))
  log(sprintf("\toutput file: %s", filenames$rds))
  log(sprintf("\tparameters file: %s", filenames$pars))
  log(sprintf("\tlog file: %s", filenames$log))

  if (file.exists(filenames$rds)) {
    bounds <- readRDS(filenames$rds)$community$bounds
  } else {
    bounds <- trait_bounds(trait_names, p)
  }

  empty <- is.null(bounds)
  if (empty) {
    bounds <- trait_bounds(trait_names, NULL)
  }

  ## Construct the assembly object
  control <- list(birth_move_tol=1, compute_viable_fitness=FALSE)
  sys0 <- plant.assembly::community(p, bounds,
    fitness_approximate_control=list(type="gp"))
  obj <- plant.assembly::assembler(sys0, control=control, filename=filenames$rds)

  ## And off we go:
  if (!empty) {
    obj <- plant.assembly::assembler_run(obj, nsteps)
  }

  invisible(obj)
}

## Utility things here:
stop_logging <- function() {
  loggr::log_info("Ending log")
  loggr::deactivate_log()
}


## Common bounds for traits; we'll use this to set maximum bounds so
## that the traits don't run out to points that make no sense.
trait_bounds <- function(trait_names, p=NULL) {
  min_bounds <- plant::bounds(lma=c(0.01, 4),
                              hmat=c(1.0, 100))
  max_bounds <- plant::bounds(lma=10^c(-5, 2),
                              hmat=c(1.0, 100))

  i <- match(trait_names, rownames(max_bounds))
  if (any(is.na(i))) {
    stop("Unknown traits: ", trait_names[is.na(i)])
  }
  if (is.null(p)) {
    bounds <- min_bounds[i, , drop=FALSE]
  } else {
    bounds <- plant::viable_fitness(max_bounds[i, , drop=FALSE], p, x=mean(min_bounds[i, , drop=FALSE]))
    bounds[, "lower"] <- pmin(bounds[, "lower"], min_bounds[i, "lower"])
    bounds[, "upper"] <- pmax(bounds[, "upper"], min_bounds[i, "upper"])
  }
  bounds
}
