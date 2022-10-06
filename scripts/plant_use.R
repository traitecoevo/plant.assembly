## This is the R version of the python code in tree_use.py
library(plant)
library(plant.assembly)

time_disturbance <- 11.3
slope <- 1.3

## These are the 'X' values:
lma <- c(0.018, 0.09485587)
## These are the 'Y' values that happen to be quite good:
birth_rate <- c(450.587114143259, 1.31521896285294)

sys <- for_python_equilibrium(time_disturbance, slope, lma, birth_rate,
                              verbose=TRUE)

w <- for_python_fitness(sys, lma)
w

lma2 <- seq_log_range(sys$bounds, 101)
ww <- for_python_fitness_approximate(sys, lma2)

plot(lma2, ww, log="x", type="l")
abline(h=0)
abline(v=sys$traits)
points(sys$traits, w)
rug(sys$fitness_approximate_points[, "lma"])

## Then, try evolution from this point:
sys2 <- for_python_evolve(time_disturbance, slope, 20,
                          sys$traits, sys$birth_rate, verbose=TRUE)

ww2 <- for_python_fitness_approximate(sys2, lma2)

## The equilibrium community vs the one from above:
plot(lma2, ww2, log="x", type="l")
lines(lma2, ww, col="grey")
abline(h=0)
abline(v=sys2$traits)
rug(sys2$fitness_approximate_points[, "lma"])

## From an empty community:
sys3 <- for_python_evolve(time_disturbance, slope, 20, verbose=TRUE)
ww3 <- for_python_fitness_approximate(sys3, lma2)

## The community assembled from scratch vs the one with the head start:
plot(lma2, ww2, log="x", type="l")
lines(lma2, ww3, col="red")
abline(h=0)
abline(v=sys2$traits)
abline(v=sys3$traits, col="red")
