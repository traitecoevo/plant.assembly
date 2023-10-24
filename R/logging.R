plant_log_info <- function(...) {
  loggr::log_info(..., package="plant", pid=Sys.getpid())
}

plant_log_debug <- function(...) {
  loggr::log_debug(..., package="plant", pid=Sys.getpid())
}

plant_log_max_fitness <- function(...) {
  plant_log_info(..., routine="max_fitness")
}

plant_log_assembler <- function(...) {
  plant_log_info(..., routine="assembler")
}

plant_log_births <- function(...) {
  plant_log_info(..., routine="births")
}

plant_log_state <- function(...) {
  plant_log_info(..., routine = "state")
}

plant_log_deaths <- function(...) {
  plant_log_info(..., routine="deaths")
}
