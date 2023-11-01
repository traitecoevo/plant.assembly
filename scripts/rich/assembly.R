devtools::load_all("..")
library(testthat)
library(plant)
plant_log_console()

p <- scm_base_parameters()
p$disturbance_mean_interval <- 7.0
sys0 <- community(p,  plant::bounds_infinite("lma"))

## This is not working with a change to the expand parameters code.
## Should possibly nuke existing cohort times to be simple.
obj_m0 <- assembler(sys0, list(birth_move_tol=0))
obj_m <- assembler_run(obj_m0, 20)
obj_m$done

ff <- lapply(obj_m$history, community_fitness_approximate)
cols <- c("black", "blue", "orange")
lma <- seq_log_range(obj_m0$community$bounds, 400)
w <- sapply(ff, function(f) f(lma))
res_lma <- sapply(obj_m$history[-1], function(x) x$traits)
res_w0 <- sapply(seq_along(res_lma), function(i) ff[[i  ]](res_lma[[i]]))
res_w1 <- sapply(seq_along(res_lma), function(i) ff[[i+1]](res_lma[[i]]))
matplot(lma, w, type="l", col=cols, lty=1, ylim=c(-1, max(w)),
        log="x")
abline(h=0, col="grey")
segments(res_lma, res_w0, res_lma, res_w1, col=cols[-1])
points(res_lma, res_w1, col=cols[-1], pch=19)

## Again, with attempting to move things:
obj_n0 <- assembler(sys0, list(birth_move_tol=1))
obj_n <- assembler_run(obj_n0, 20)

ff <- lapply(obj_n$history, community_fitness_approximate)
w <- sapply(ff, function(f) f(lma))
res_lma <- sapply(obj_n$history[-1], function(x) x$traits)
res_w0 <- sapply(seq_along(res_lma), function(i) ff[[i  ]](res_lma[[i]]))
res_w1 <- sapply(seq_along(res_lma), function(i) ff[[i+1]](res_lma[[i]]))
matplot(lma, w, type="l", col=cols, lty=1, ylim=c(-1, max(w)),
        log="x")
abline(h=0, col="grey")
segments(res_lma, res_w0, res_lma, res_w1, col=cols[-1])
points(res_lma, res_w1, col=cols[-1], pch=19)

## And stochastic:
set.seed(1)
p <- scm_base_parameters()
p$disturbance_mean_interval <- 7.0
sys0 <- community(p,  plant::bounds_infinite("lma"))
max_bounds <- bounds(lma=c(0.01, 10))
vcv <- mutational_vcv_proportion(max_bounds, 0.001)
obj_s0 <- assembler(sys0,
                    list(birth_type="stochastic",
                         run_type="single",
                         vcv=vcv))

obj_s <- assembler_run(obj_s0, 30)

tmp <- community_prepare_approximate_fitness(obj_s$community)
ff <- community_fitness_approximate(tmp)
w <- ff(lma)
plot(lma, w, log="x", ylim=c(-1, max(w)), type="l")
points(obj_s$community$traits, rep(0, length(obj_s$community$traits)))
