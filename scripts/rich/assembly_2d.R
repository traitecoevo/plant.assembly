devtools::load_all("..")
library(testthat)
library(plant)
plant_log_console()

assembler_parameters <- function(time_disturbance, B4=2.0, B5=2.0) {
  p <- scm_base_parameters()
  p$disturbance_mean_interval <- time_disturbance
  p$hyperpar <- make_FF16_hyperpar(B4=B4, B5=B5)
  p
}

max_bounds <- bounds(lma=10^c(-4, 2), rho=c(1, 1e4))
p <- assembler_parameters(10)
# bounds <- plant::viable_fitness(max_bounds, p)

## Jump straight to the end:
bounds <- bounds(lma=c(0.0703621360996001, 0.938239530095711),
                 rho=c(54.9432041575792, 608))

sys0 <- community(p, bounds)
obj_m0 <- assembler(sys0, list(birth_move_tol=1, compute_viable_fitness=FALSE))
if (file.exists("obj_m.rds")) {
  obj_m <- readRDS("obj_m.rds")
  saveRDS(obj_m, "obj_m.rds")
} else {
  obj_m <- assembler_run(obj_m0, 20)
}

lscape <- parallel::mclapply(obj_m$history[-1], local_landscape)

zmin <- -.5
zmax <- 5
xlim <- range(sapply(lscape, function(x) range(x$x)))
ylim <- range(sapply(lscape, function(x) range(x$y)))
zlim <- c(-zmax, zmax)

pdf("assembly_2d.pdf")
for (i in lscape) {
  plot(i, xlim=xlim, ylim=ylim, zlim=zlim, zmin=zmin, zmax=zmax)
}
dev.off()
