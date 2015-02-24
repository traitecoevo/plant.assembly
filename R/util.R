copy_attributes <- function(from, to, exclude=character(0)) {
  at_from <- attributes(from)
  for (i in setdiff(names(at_from), c("dim", "dimnames"))) {
    attr(to, i) <- at_from[[i]]
  }
  to
}

## TODO: Get some testing in here so I can drop the asserts and merge
## with closest2
##
## TODO:
closest <- function(pt, A, r=NULL) {
  ## stopifnot(ncol(A) == length(pt))
  if (!is.null(r)) {
    ## stopifnot(is.matrix(r) &&
    ##                         ncol(r) == 2L &&
    ##                         nrow(r) == ncol(A))
    dr <- r[,2] - r[,1]
    A <- t((t(A) - r[,1]) / dr)
    pt <- (pt - r[,1]) / dr
    ## TODO:?
    ## A <- rescale(A, r)
    ## pt <- rescale(pt, r)
  }
  d <- colSums((t(A) - drop(pt))^2)
  i <- which.min(d)
  attr(i, "distance") <- sqrt(d[i])
  i
}

closest_log <- function(pt, A, r=NULL) {
  closest(log(pt), log(A), if (is.null(r)) NULL else log(r))
}

## Possible rewrite?
## closest2 <- function(pt, A) {
##   d <- colSums((t(A) - drop(pt))^2)
##   i <- which.min(d)
##   attr(i, "distance") <- sqrt(d[i])
##   i
## }

## Same rescaling issue here as for determining distance to points in
## the maximisation.
rescale <- function(A, r=NULL) {
  if (is.null(r)) {
    r <- t(apply(A, 2, range))
    colnames(r) <- c("lower", "upper")
  }
  stopifnot(is.matrix(r) &&
                          ncol(r) == 2L &&
                          nrow(r) == ncol(A))
  dr <- r[,2] - r[,1]
  t((t(A) - r[,1]) / dr)
}

unrescale <- function(A, r) {
  stopifnot(is.matrix(r) &&
            ncol(r) == 2L &&
            nrow(r) == ncol(A))
  dr <- r[,2] - r[,1]
  t(t(A) * dr + r[,1])
}

maximize_logspace <- function(f, x, bounds, tol) {
  control <- list(tol=tol, maximize=TRUE)
  log_bounds <- log(bounds)
  fit <- dfoptim::nmkb(log(x), function(x) f(exp(x)),
                       lower=log_bounds[,1], upper=log_bounds[,2],
                       control=list(tol=tol, maximize=TRUE))
  fit$par <- exp(fit$par)
  fit
}

has_attr <- function(x, which) {
  !is.null(attr(x, which, exact=TRUE))
}
