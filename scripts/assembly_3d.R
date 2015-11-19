devtools::load_all("..")
library(testthat)
library(plant)
plant_log_console()

## NOTE: This is not using the "hybrid" equilibrium solver, but I
## think that it should be.
assembler_parameters <- function(time_disturbance, B4=2.0, B5=2.0) {
  p <- scm_base_parameters()
  p$disturbance_mean_interval <- time_disturbance
  p$hyperpar <- make_FF16_hyperpar(B4=B4, B5=B5)
  p$control$equilibrium_solver_name <- "hybrid"
  p
}

## Times are 2, 5, 10, 20, 40
max_bounds <- rbind(lma=10^c(-4, 2),
                    rho=c(1, 1e4),
                    hmat=c(1.0, 100))

f <- function(t) {
  p <- assembler_parameters(t)
  ## bounds <- viable_fitness(max_bounds, p)

  ## Jump straight to the end:
  bounds <- bounds(lma=c(0.0273839008086021, 1.33352391701713),
                   rho=c(3.16230123030217, 9999.59628838755),
                   hmat=c(1, 16.5958691))

  filename <- sprintf("3d_pilot_%d.rds", p$disturbance_mean_interval)
  control <- list(birth_move_tol=1, compute_viable_fitness=FALSE)

  sys0 <- community(p, bounds)
  obj_m0 <- assembler(sys0, control=control, filename=filename)
  assembler_run(obj_m0, 20)
}

times <- c(2, 4, 8, 16, 32, 64, 128)
res <- parallel::mclapply(times, f, mc.preschedule=FALSE, mc.cores=3L)

if (file.exists("obj_m3.rds")) {
  obj_m3 <- readRDS("obj_m3.rds")
} else {
  obj_m3 <- assembler_run(obj_m0, 20)
  saveRDS(obj_m3, "obj_m3.rds")
}

obj <- readRDS(filename)

## TODO: looks like we're saving the wrong thing:
lscape <- parallel::mclapply(obj[-1], local_landscape_matrix)

## Compute all the limits
lim <- lapply(lscape, function(x) do.call("rbind", diag(x)))
lim <- t(apply(do.call("cbind", lim), 1, range))

zmin <- -2
zmax <- 5
zlim <- c(-zmax, zmax)

log <- c(TRUE, TRUE, FALSE)

pdf(sprintf("3d_pilot_%02d.pdf", p$disturbance_mean_interval))
for (m in lscape) {
  plot(m, lim, zlim=zlim, zmin=zmin, zmax=zmax, log=log)
}
dev.off()
