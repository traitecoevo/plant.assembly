##' Generate a local landscape from a history with a
##' \code{fitness_approximate_slopes} element.
##' @title Generate local landscape
##' @param sys Single community with element
##'   \code{fitness_approximate_slopes}
##' @param xlim,ylim Optional x and y limits (otherwise resident range
##'   will be expanded by \code{scal})
##' @param scal Scaling parameter when creating \code{xlim} and
##'   \code{ylim}.
##' @param n Number of points to compute fitness for (in each of x and
##'   y, so \code{n^2} points total are used)
##' @param combine Leave this be.
##' @param traits Which traits to compute the matrix for.  Not
##'   required with two traits, but for more than two traits, you must
##'   indicate which pair of traits to compute the landscape for.
##' @export
local_landscape <- function(sys, xlim=NULL, ylim=NULL, scal=.25,
                            n=101, combine="closest", traits=NULL) {
  X <- sys$traits
  if (ncol(X) < 2L) {
    stop("This requires at least two dimensions")
  }

  if (is.null(traits)) {
    traits <- c(1L, 2L)
  } else if (is.character(traits)) {
    traits <- match(traits, colnames(X))
    if (any(is.na(traits))) {
      stop("Unknown traits")
    }
  } else if (any(traits > ncol(X))) {
    stop("Traits out of bounds")
  }

  if (length(traits) != 2L) {
    stop("Expected two traits")
  }

  scal <- 1 + scal
  if (is.null(xlim)) {
    xlim <- range(X[, traits[[1]]]) * c(1/scal, scal)
  }
  if (is.null(ylim)) {
    ylim <- range(X[, traits[[2]]]) * c(1/scal, scal)
  }

  f <- make_approximate_fitness_slopes(sys, combine, traits)
  x <- seq_log_range(xlim, n)
  y <- seq_log_range(ylim, n)
  xy <- as.matrix(expand.grid(x, y))
  colnames(xy) <- sys$trait_names[traits]
  z <- matrix(apply(xy, 1, f), n, n)
  ret <- list(x=x, y=y, xy=xy, z=z, resident=X[, traits, drop=FALSE])
  class(ret) <- "local_landscape"
  ret
}

local_landscape_matrix <- function(sys, ...) {
  n <- length(sys$trait_names)
  m <- vector("list", n^2)
  dim(m) <- c(n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(i - 1L)) {
      lscape <- local_landscape(sys, traits=c(j, i), ...)
      m[[j, i]] <- lscape
      m[[i, j]] <- list(x=lscape$y, y=lscape$x,
                        xy=lscape$xy[,2:1,drop=FALSE],
                        z=t(lscape$z),
                        resident=lscape$resident[,2:1,drop=FALSE])
      class(m[[i, j]]) <- class(m[[j, i]])
    }
  }

  for (i in seq_len(n)) {
    if (i == 1L) {
      m[[i, i]] <- range(m[[i + 1L, i]]$y)
    } else {
      m[[i, i]] <- range(m[[i, i - 1L]]$x)
    }
  }
  rownames(m) <- colnames(m) <- sys$trait_names
  class(m) <- c("local_landscape_matrix", "matrix")
  m
}

## Useful values of combine:
##   max: take the best fitness
##   "closest": take the closest point
##   FALSE: don't combine and instead return all the points
make_approximate_fitness_slopes <- function(sys, combine=max, i=NULL) {
  xx <- log(sys$traits)

  if (is.null(i)) {
    funcs <- lapply(sys$fitness_approximate_slopes, grader::taylor2)
  } else {
    trim <- function(obj) {
      list(x=obj$x[i], fx=obj$fx, gr=obj$gr[i], H=obj$H[i, i, drop=FALSE])
    }
    funcs <- lapply(sys$fitness_approximate_slopes,
                    function(x) grader::taylor2(trim(x)))
    xx <- xx[, i, drop=FALSE]
  }

  if (identical(combine, "closest")) {
    function(x) {
      lx <- log(x)
      funcs[[closest(lx, xx)]](lx)
    }
  } else {
    if (identical(combine, FALSE)) {
      combine <- t
    } else {
      combine <- match.fun(combine)
    }
    function(x) {
      lx <- log(x)
      combine(sapply(funcs, function(f) f(lx)))
    }
  }
}

local_landscape_zlim <- function(x, zmin, zmax) {
  c(-1, 1) * max(abs(pmin(pmax(x$z, zmin), zmax)))
}

##' @export
plot.local_landscape <- function(x, xlim=NULL, ylim=NULL, zlim=NULL,
                                 zmin=-.5, zmax=5, ...) {
  if (is.null(xlim)) {
    xlim <- range(x$x)
  }
  if (is.null(ylim)) {
    ylim <- range(x$y)
  }
  x$z[x$z < zmin] <- NA
  x$z[x$z > zmax] <- zmax
  if (is.null(zlim)) {
    zlim <- c(-1, 1) * max(abs(x$z), na.rm=TRUE)
  }

  cols <- rev(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
                "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd",
                "#5e4fa2"))

  nms <- colnames(x$resident)

  image(x, log="xy", xlim=xlim, ylim=ylim, col=cols, zlim=zlim,
        las=1, xlab=nms[[1]], ylab=nms[[2]], ...)
  contour(x, levels=0, add=TRUE, labels="")
  points(x$resident, pch=19)
}

##' @export
plot.local_landscape_matrix <- function(m, lim, ..., gap=1, log=TRUE) {
  n <- ncol(m)
  labels <- colnames(m)

  opar <- par(mfrow=c(n, n), mar=rep.int(gap/2, 4), oma=rep(4, 4))
  on.exit(par(opar))

  xl <- yl <- logical(n)
  if (is.numeric(log)) {
    xl[log] <- yl[log] <- TRUE
  } else if (is.logical(log)) {
    xl[] <- yl[] <- rep(log, length.out=n)
  } else {
    xl[] <- grepl("x", log)
    yl[] <- grepl("y", log)
  }

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      l <- paste0(ifelse(xl[j], "x", ""), ifelse(yl[i], "y", ""))
      plot(NA, type="n", xlim=lim[j,], ylim=lim[i,],
           xlab="", ylab="", axes=FALSE, log=l)
      if (i == 1 && j %% 2L == 0) {
        axis(3L)
      } else if (i == n && j %% 2L != 0) {
        axis(1L)
      }
      if (j == 1 && i %% 2L == 0) {
        axis(2L)
      } else if (j == n && i %% 2L != 0) {
        axis(4L)
      }
      if (i == j) {
        par(usr=c(0, 1, 0, 1))
        xlp <- if (xl[i]) 10^0.5 else 0.5
        ylp <- if (yl[j]) 10^0.5 else 0.5
        cex_labels <- max(0.8, min(2, 0.9 / max(strwidth(labels, "user"))))
        text(xlp, ylp, labels[i], cex=cex_labels)
      } else {
        plot(m[[j, i]], add=TRUE, ...)
      }
      box()
    }
  }
}
