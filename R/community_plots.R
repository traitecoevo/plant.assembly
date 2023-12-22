#' Plot fitness landscape
#'
#' @param community XXXX
#' @param label XXX 
#' @param xlim XXX
#' @param ylim XXX
#'
#' @return community_landscape returns one plot, history_landscape returns plots length of timeseries
#' @export
community_plot_fitness_landscape <- function(community, label = NA, xlim = community$bounds, ylim = c(-20, 20)){

  trait_names <- community$trait_names
  landscape <- community$fitness_points
  landscape[["x"]] = landscape[[trait_names]]

  data_residents <-
    tibble(
      x = community$traits[,trait_names], 
      fitness = approx(landscape[[trait_names]], landscape$fitness, community$traits[, trait_names])$y
    )
  
  p <- 
    landscape %>%
    ggplot(aes(x, fitness, col=batch_nr)) +
    geom_point() +
    geom_line() +
    geom_point(data = data_residents, col = "red") +
    geom_abline(intercept = 0, slope = 0, linetype = "dashed") +
    scale_x_log10(limits = xlim) +
    ylab("Fitness") +
    #scale_y_continuous(limits = ylim) +
    #ylim(ylim) +
    theme_classic() + 
    theme(text = element_text(size = 16),
          legend.position = "none")
  
  if(!is.na(label)){
    p +
      geom_text(aes(x = xlim[1], y = ylim[2], label = label, vjust = "inward", hjust = "inward"), size = 5) -> p
  }
  
  p
}



#' Returns birth rates for residents in community
#'
#' @param tidy_community Community or history object from tidied assembly
#' @param ... 
#'
#' @return Returns one plot if community of, if history, plots length of history timeseries
#' @export
#'
plot_community <- function(tidy_community, ...){
  if(is_tibble(tidy_community)){plot_community_1d(tidy_community, ...) -> p}
  if(!is_tibble(tidy_community)){purrr::imap(tidy_community, ~plot_community_1d(tidy_community = .x, step = .y), ...)-> p}
  
  invisible(p)
  
}

plot_community_1d <- function(tidy_community, step = NA, xlim = c(0.01, 1), ylim = c(1e-4, 5)){
  tidy_community %>%
    select(-births, -invader) %>%
    names() -> traits
  tidy_community %>%
    ggplot(aes_string(x = traits[1], y = "births")) + 
    geom_point(aes(colour = invader), size = 2) +
    xlab(traits[1]) +
    ylab("Birth rate") +
    theme_classic() + 
    theme(text = element_text(size = 16),
          legend.position = "none") +
    scale_x_log10(limits = xlim) +
    scale_y_log10(limits = ylim) -> p 
  
  if(!is.na(step)){
    p +
      geom_text(aes(x = xlim[1], y = ylim[2], label = paste0("Step = ",step), vjust = "inward", hjust = "inward"), size = 5) -> p
  }
  
  return(p)
}

#' Plot pair-wise trait combinations of resident community 
#'
#' @param tidy_community Community or history object from tidied assembly with two traits
#' @param ... 
#'
#' @return Returns one plot if community of, if history, plots length of history timeseries
#' @export
#'
plot_community_2d <- function(tidy_community, ...){
  if(is_tibble(tidy_community)){
    
    ylim_max <- max(tidy_community$hmat)
    ylim_min <- min(tidy_community$hmat)
    ylim <- c(ylim_min, ylim_max)
    
    xlim_max <- max(tidy_community$lma)
    xlim_min <- min(tidy_community$lma)
    xlim <- c(xlim_min, xlim_max)
    
    plot_community_2d_internal(tidy_community, xlim = xlim, ylim = ylim, ...) -> p
  }
  if(!is_tibble(tidy_community)){
    ylim_max <- max(map_dbl(tidy_community, ~max(.x$hmat)))
    ylim_min <- min(map_dbl(tidy_community, ~min(.x$hmat)))
    ylim <- c(ylim_min, ylim_max)
    
    xlim_max <- max(map_dbl(tidy_community, ~max(.x$lma)))
    xlim_min <- min(map_dbl(tidy_community, ~min(.x$lma)))
    xlim <- c(xlim_min, xlim_max)
    
    purrr::imap(tidy_community, ~plot_community_2d_internal(tidy_community = .x, step = .y, xlim = xlim, ylim = ylim), ...)-> p
  }
  
  invisible(p)
  
}

plot_community_2d_internal <- function(tidy_community, step = NA, xlim = c(0.01, 1), ylim = c(1e-4, 5)){
  tidy_community %>%
    select(-births, -invader) %>%
    names() -> traits
  tidy_community %>%
    ggplot(aes_string(x = traits[1], y = traits[2])) + 
    geom_point(aes(colour = births), size = 2) +
    theme_classic() + 
    theme(text = element_text(size = 16),
          legend.position = "none") +
    scale_x_log10(limits = xlim) +
    scale_y_log10(limits = ylim) -> p 
  
  if(!is.na(step)){
    p +
      geom_text(aes(x = xlim[1], y = ylim[2], label = paste0("Step = ",step), vjust = "inward", hjust = "inward"), size = 5) -> p
  }
  
  return(p)
}
