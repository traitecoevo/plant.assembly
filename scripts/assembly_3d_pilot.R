## Running this takes a long time; around 2 weeks I think.  So it'll
## move from this directory and into successional diversity.
##
## TODO: The assembler would ideally be self-restarting which means
## control and filename (etc) needs to be stored within the returned
## object.
devtools::load_all("../../plant")
library(testthat)
library(plant)
library(plant.ml)
library(loggr)
library(parallel)
plant_log_console()

path <- "output"
dir.create(path, FALSE, TRUE)

## NOTE: This is not using the "hybrid" equilibrium solver, but I
## think that it should be.
assembler_parameters <- function(time_disturbance, B4=2.0, B5=2.0) {
  p <- ebt_base_parameters()
  p$disturbance_mean_interval <- time_disturbance
  p$hyperpar <- make_FFW16_hyperpar(B4=B4, B5=B5)
  p$control$equilibrium_solver_name <- "hybrid"
  p
}

filename_simulation <- function(t) {
  file.path(path, sprintf("3d_pilot_%d.rds", t))
}

run_simulation <- function(t, steps) {
  ## These are values that need to be included:
  min_bounds <- bounds(lma=c(0.01, 4),
                       rho=c(10, 1000),
                       hmat=c(1.0, 100))
  max_bounds <- bounds(lma=10^c(-4, 2),
                       rho=c(1, 1e4),
                       hmat=c(1.0, 100))
  p <- assembler_parameters(t)
  filename <- filename_simulation(t)
  control <- list(birth_move_tol=1, compute_viable_fitness=FALSE)

  if (file.exists(filename)) {
    ## bounds <- last(readRDS(filename))$bounds
    bounds <- max_bounds
  } else {
    ## could easily move the viable fitness bit within compute_viable_bounds
    bounds <- viable_fitness(max_bounds, p)
    bounds[, "lower"] <- pmin(bounds[, "lower"], min_bounds[, "lower"])
    bounds[, "upper"] <- pmax(bounds[, "upper"], min_bounds[, "upper"])
  }

  sys0 <- community(p, bounds)
  obj_m0 <- assembler(sys0, control=control, filename=filename)
  assembler_run(obj_m0, steps)
}

times <- c(2, 4, 8, 16, 32, 64)
if (FALSE) {
  res <- parallel::mclapply(rev(times), run_simulation, 60,
                            mc.preschedule=FALSE, mc.cores=3L)
}

dat <- lapply(filename_simulation(times), readRDS)
dat_last <- lapply(dat, last)

## Computing "true" landscapes takes quite a while; expect this to
## take a while.
lscape_approx <- lapply(dat_last, local_landscape_matrix)
if (FALSE) {
  lscape_real <- mclapply(dat_last, local_landscape_matrix,
                          control=list(method="real", n_real=50L),
                          mc.preschedule=FALSE)
  saveRDS(lscape_real, file.path(path, "lscape_real.rds"))
}
lscape_real <- readRDS(file.path(path, "lscape_real.rds"))

pdf(file.path(path, "landscape_approx.pdf"))
for (x in lscape_approx) {
  plot(x)
}
dev.off()

pdf(file.path(path, "landscape_real.pdf"))
for (x in lscape_real) {
  plot(x)
}
dev.off()
