devtools::load_all("..")
library(testthat)
library(plant)
plant_log_console()

assembler_parameters <- function(time_disturbance, B4=2.0, B5=2.0) {
  p <- ebt_base_parameters()
  p$disturbance_mean_interval <- time_disturbance
  p$hyperpar <- make_FFW16_hyperpar(B4=B4, B5=B5)
  p
}

max_bounds <- bounds(lma=10^c(-4, 2), rho=c(1, 1e4))
p <- assembler_parameters(10)
# bounds <- viable_fitness(max_bounds, p)

## Jump straight to the end:
bounds <- bounds(lma=c(0.0703621360996001, 0.938239530095711),
                 rho=c(54.9432041575792, 608))

sys0 <- community(p, bounds)
obj_m0 <- assembler(sys0, list(birth_move_tol=1, compute_viable_fitness=FALSE))
obj_m <- assembler_run(obj_m0, 20)
