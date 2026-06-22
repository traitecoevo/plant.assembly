# Internal finite-difference gradient helpers.
#
# Vendored from richfitz/grader (gradient_points / gradient_extrapolate) so the
# package does not depend on that tiny, non-CRAN package. Used by
# community_selection_gradient() in solve_attractors.R.
#
# gradient_points() builds the set of evaluation points (a sequence of
# successively halved step sizes per dimension) needed for a central-difference
# estimate with optional Richardson extrapolation. gradient_extrapolate() turns
# the fitness values at those points into the gradient.

gradient_points <- function(x, eps = 1e-04, d = 1e-04, r = 4,
                            zero_tol = sqrt(.Machine$double.eps / 7e-07)) {
  v <- 2
  n <- length(x)
  h <- abs(d * x) + eps * (abs(x) < zero_tol)
  pts <- vector("list", r * n)
  dim(pts) <- c(r, n)
  dx <- matrix(NA, r, n)
  j <- seq_len(n)
  for (k in seq_len(r)) {
    for (i in seq_len(n)) {
      dx_i <- h * (i == j)
      pts[[k, i]] <- rbind(x + dx_i, x - dx_i)
    }
    dx[k, ] <- 2 * h
    h <- h / v
  }
  ret <- do.call("rbind", pts)
  attr(ret, "dim_y") <- c(2L, r * n)
  attr(ret, "dx") <- dx
  attr(ret, "n") <- n
  attr(ret, "r") <- r
  class(ret) <- "gradient_points"
  ret
}

gradient_extrapolate <- function(y, pts) {
  dx <- attr(pts, "dx")
  r <- attr(pts, "r")
  n <- attr(pts, "n")
  a <- (y[1, ] - y[2, ]) / dx
  for (m in seq_len(r - 1L)) {
    four_m <- 4^m
    a_next <- matrix(NA, r - m, n)
    for (i in seq_len(nrow(a_next))) {
      a_next[i, ] <- (a[i + 1L, ] * four_m - a[i, ]) / (four_m - 1)
    }
    a <- a_next
  }
  drop(a)
}
