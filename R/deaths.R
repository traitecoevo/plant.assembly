assembler_deaths <- function(obj) {
  obj$community <- community_deaths(obj$community, obj$control)
  obj
}

## Really simple deaths.
community_deaths <- function(community, control) {
  eps <- control$dead_seed_rain
  check_inviable <- control$check_inviable
  to_drop <- community$seed_rain < eps
  if (any(to_drop)) {
    plant_log_deaths(paste0("Dropping species ",
                            paste(which(to_drop), collapse=", ")))
    community <- community_drop(community, to_drop)
  }
  if (check_inviable && length(community) > 1L) {
    res <- plant::check_inviable(community_parameters(community))
    community$seed_rain <- as.numeric(res)
    to_drop2 <- attr(res, "drop")
    if (any(to_drop2)) {
      plant_log_deaths(paste0("Dropping inviable species ",
                              paste(which(to_drop2), collapse=", ")))
      community <- community_drop(community, to_drop2)
    }
  }
  community
}
