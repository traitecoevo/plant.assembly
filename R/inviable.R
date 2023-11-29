##' Check low-abundance strategies for viability.
##'
##' @title Check low-abundance strategies for viability
##' @param p A Parameters object
##' @param ctrl Control object
##' @export
check_inviable <- function(p, ctrl) {
  ## eps_test: *Relative* value to use for determining what
  ## "low abundance" means.  Species that have a offspring arrival of less than
  ## `eps_test * max(p$birth_rate)` will be tested.  By default
  ##  this is 1 100th of the maximum offspring arrival.
  ## TODO: don't do anything if we don't have at least 2 species?
  eps <- ctrl$equilibrium_extinct_birth_rate
  ## TODO: This was ctrl$equilibrium_inviable_test, but I think
  ## that birth offspring arrival actually makes more sense?  It's fractional
  ## though so who knows.
  eps_test <- 1e-2
  birth_rate <- sapply(p$strategies, function(s) s$birth_rate_y, simplify = TRUE)
  ## NOTE: We don't actually run to equilibrium here; this is just
  ## because it's a useful way of doing incoming -> outgoing offspring
  ## rain.
  runner <- make_equilibrium_runner(p, ctrl =ctrl)
  offspring_production <- runner(birth_rate)
  
  test <- which(offspring_production < birth_rate &
                  birth_rate < max(offspring_production) * eps_test)
  test <- test[order(offspring_production[test])]
  
  drop <- logical(length(offspring_production))
  
  for (i in test) {
    plant_log_inviable(paste("Testing species", i),
                       stage="testing", species=i)
    x <- offspring_production
    x[i] <- eps
    res <- runner(x)
    if (res[[i]] < eps) {
      plant_log_inviable(paste("Removing species", i),
                         stage="removing", species=i)
      drop[[i]] <- TRUE
      res[[i]] <- 0.0
      offspring_production <- res
    }
  }
  
  ## It's possible that things slip through and get driven extinct by
  ## the time that they reach here.
  drop <- drop | offspring_production < eps
  
  attr(offspring_production, "drop") <- drop
  offspring_production
}
