## The fitness code.  This is really important and hard to get right.
## There are a few things here:
devtools::load_all("..")
library(testthat)
library(plant)
plant_log_console()

p <- ebt_base_parameters()
sys0 <- community(p, bounds_infinite("lma"))
sys0 <- community_viable_bounds(sys0)

sys0_gp <- community(p, sys0$bounds,
                     fitness_approximate_control=list(type="gp"))

sys0 <- community_prepare_approximate_fitness(sys0)
sys0_gp <- community_prepare_approximate_fitness(sys0_gp)

f <- community_fitness_approximate(sys0)
f_gp <- community_fitness_approximate(sys0_gp)

xx <- seq_log_range(sys0$bounds, 500)

plot(sys0$fitness_approximate_points, log="x")
points(sys0_gp$fitness_approximate_points, col="blue")

lines(xx, f(xx))
lines(xx, f_gp(xx), col="blue")

plot(xx, f(xx) - f_gp(xx), log="x", type="l")

## Then, using this to create a births function...
g <- function(sys) {
  community_new_types_maximum_fitness(sys, assembler_control(NULL))
}

res <- g(sys0)
res_gp <- g(sys0_gp)

plot(xx, f(xx), log="x", type="l")
points(res, attr(res, "fitness"), pch=19)
points(res_gp, attr(res_gp, "fitness"), col="red", cex=2)
