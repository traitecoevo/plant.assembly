##' Helper function for making bounds
##' @title Trait bounds
##' @param ... Named list, each element of which is a 2-element
##' numeric vector of lower and upper bounds.
##' @export
##' @examples
##' bounds(lma=c(0.01, 10))
##' bounds(lma=c(0.01, 10), rho=c(1, 1000))
bounds <- function(...) {
  x <- list(...)
  if (length(x) == 0) {
    stop("Need at least one argument")
  }
  if (!all(vapply(x, length, integer(1)) == 2)) {
    stop("All entries must be length 2")
  }
  if (is.null(names(x)) || any(names(x) == "")) {
    stop("All elements must be named")
  }
  ret <- rbind_list(x)
  colnames(ret) <- c("lower", "upper")
  ret
}

##' @param bounds A set of bounds
##' @param finite Logical indicating if bounds must be finite
##' @rdname bounds
##' @export
check_bounds <- function(bounds, finite=FALSE) {
  if (!is.matrix(bounds)) {
    stop("bounds must be a matrix")
  }
  if (ncol(bounds) != 2) {
    stop("bounds must have two columns")
  }
  if (is.null(rownames(bounds))) {
    stop("bounds must have rownames")
  }
  if (finite && any(!is.finite(bounds))) {
    stop("bounds must be finite")
  }
  colnames(bounds) <- c("lower", "upper")
  invisible(bounds)
}

##' @param trait_names Character vector of trait names
##' @rdname bounds
##' @export
bounds_infinite <- function(trait_names) {
  n <- length(trait_names)
  b <- cbind(lower=rep(-Inf, n), upper=rep(Inf, n))
  rownames(b) <- trait_names
  b
}

##' @export
##' @rdname bounds
##' @param x a point to detect if it lies within bounds
check_point <- function(x, bounds) {
  if (is.matrix(x)) {
    if (ncol(x) != nrow(bounds)) {
      stop("Invalid size x")
    }
  } else {
    if (length(x) != nrow(bounds)) {
      stop("Invalid size x")
    }
    x <- rbind(x, deparse.level=0)
  }
  if (is.null(names(x))) {
    colnames(x) <- rownames(bounds)
  } else if (names(x) != rownames(bounds)) {
    stop("Incorrect names on x")
  }
  tx <- t(x)
  if (any(tx < bounds[, "lower"] | tx > bounds[, "upper"])) {
    stop("Value does not lie within bounds")
  }
  invisible(x)
}
