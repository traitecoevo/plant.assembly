##' Generate a local landscape from a history with a
##' \code{fitness_approximate_slopes} element.
##' @title Generate local landscape
##' @param sys Single community with element
##'   \code{fitness_approximate_slopes}
##' @param xlim,ylim Optional x and y limits (otherwise resident range
##'   will be expanded by \code{control$scal})
##' @param control Set of control parameters, passed to
##'   \code{local_landscape_control}
##' @param traits Which traits to compute the matrix for.  Not
##'   required with two traits, but for more than two traits, you must
##'   indicate which pair of traits to compute the landscape for.
##' @export
local_landscape <- function(sys, xlim=NULL, ylim=NULL, control=NULL,
                            traits=NULL) {
  control <- local_landscape_control(control)
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

  scal <- 1 + control$scal
  if (is.null(xlim)) {
    xlim <- range(X[, traits[[1]]]) * c(1/scal, scal)
  }
  if (is.null(ylim)) {
    ylim <- range(X[, traits[[2]]]) * c(1/scal, scal)
  }

  n <- control$n
  x <- seq_log_range(xlim, n)
  y <- seq_log_range(ylim, n)
  xy <- as.matrix(expand.grid(x, y))
  colnames(xy) <- sys$trait_names[traits]
  resident <- X[, traits, drop=FALSE]

  ## Roll method, zlim, combine etc together into one control object.
  if (control$method == "slopes") {
    z <- local_landscape_slopes(sys, xy, n, control)
  } else if (control$method == "real") {
    z <- local_landscape_true(sys, xy, n, control)
  } else {
    stop("Unknown method: ", control$method)
  }
  ret <- list(x=x, y=y, xy=xy, z=z, resident=resident)
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

  nms <- colnames(x$resident)

  image(x, log="xy", xlim=xlim, ylim=ylim, col=cols_landscape(),
        zlim=zlim, las=1, xlab=nms[[1]], ylab=nms[[2]], ...)
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

local_landscape_slopes <- function(sys, xy, n, control) {
  traits <- colnames(xy)
  f <- make_approximate_fitness_slopes(sys, "closest", traits)
  matrix(apply(xy, 1, f), n, n)
}

## It might make sense to compute this the way I am currently doing it
## and then resample using something like akima::interp or
## fields::interp.surface; that would allow a lower true data
## computing overhead while giving nice looking plots.
##' @importFrom progress progress_bar
local_landscape_true <- function(sys, xy, n, control) {
  zmin <- control$zmin
  n_batch <- control$n_batch
  traits <- colnames(xy)
  browser()
  ## Convert dimensions onto 0..1:
  lim <- log(t(apply(xy, 2, range)))
  resident <- rescale(log(sys$traits[, traits, drop=FALSE]), lim)
  xy_orig <- xy
  xy <- rescale(log(xy), lim)

  ## Identify the resident that the points are closest to (this
  ## actually takes a little while and would be much more efficient
  ## using a triangulation approach), which I think I actually have
  ## implemented somewhere.
  k <- apply(xy, 1, closest, resident)

  d <- sqrt(rowSums((xy - resident[k, , drop=FALSE])^2))

  ## Various indices we'll track:
  ord <- order(d)
  excl <- done <- integer(0)
  w <- matrix(NA_real_, n, n)

  make_traits <- function(i) {
    ret <- sys$traits[k[i], , drop=FALSE]
    ret[, traits] <- xy_orig[i, , drop=FALSE]
    ret
  }

  drop_ray <- function(i) {
    x0 <- resident[k[i], ]
    x1 <- xy[i, ]
    pos <- intersect(ord, which(k == k[i]))

    if (length(pos) > 0L) {
      theta <- atan2(x1[2] - x0[2],          x1[1] - x0[1])
      alpha <- atan2(sqrt(2) / (n - 1L) / 2, sqrt(sum((x1 - x0)^2))) * 1.01
      angle <- atan2(xy[pos, 2] - x0[2],     xy[pos, 1] - x0[1])

      if (theta + alpha > pi) {
        angle[angle < - pi / 2] <- angle[angle < - pi / 2] + 2 * pi
      } else if (theta - alpha < -pi) {
        angle[angle > pi / 2]   <- angle[angle >   pi / 2] - 2 * pi
      }

      pos <- pos[angle > (theta - alpha) & angle < (theta + alpha)]
    }
    pos
  }

  p <- progress::progress_bar$new(total=length(ord))$tick
  p(0)

  f <- community_make_fitness(sys)
  while (length(ord) > 0L) {
    idx <- ord[seq_len(min(n_batch, length(ord)))]

    w[idx] <- wi <- f(make_traits(idx))
    done <- c(done, idx)
    ord <- setdiff(ord, idx)
    p(length(idx))

    if (any(wi < zmin)) {
      for (idx2 in idx[wi < zmin]) {
        drop <- drop_ray(idx2)
        ord <- setdiff(ord, drop)
        p(length(drop))
        excl <- c(excl, drop)
      }
    }
  }

  w
}

cols_landscape <- function() {
  rev(c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08b",
        "#ffffbf", "#e6f598", "#abdda4", "#66c2a5", "#3288bd",
        "#5e4fa2"))
}

##' @export
##' @rdname local_landscape
local_landscape_control <- function(control) {
  defaults <- list(scal=0.25,   # limit expansion factor
                   n=101L,      # grid size
                   method="slopes",
                   ## For _real:
                   n_batch=20L, # batch size
                   zmin=-1)     # z axis cut off
  if (length(control) > 0L && is.null(names(control))) {
    stop("control must be named")
  }
  extra <- setdiff(names(control), names(defaults))
  if (length(extra) > 0L) {
    stop("Unknown entries in control: ", paste(extra, collapse=", "))
  }
  ret <- modifyList(defaults, as.list(control))
  if (!(ret$method %in% c("slopes", "real"))) {
    stop("Unknown method: ", ret$method)
  }
  ret
}
