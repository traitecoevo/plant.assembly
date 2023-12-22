#' Tidies output from assembler_run.
#'
#' @param assembly 
#'
#' @return Returns the birth rate of the final community as well as the historical timeseries of birth rates for residents. If fitness landscapes are available, provides those as well.
#' @export
#'
tidy_assembly <- function(assembly){
  
  history <- tidy_history_community(assembly$history)
  community <- tidy_community(assembly$community)
  history_landscape <- NA
  community_landscape <- NA
  if(!is.null(assembly$community$fitness_points)){
    history_landscape <- tidy_landscape(assembly$history)
    community_landscape <- tidy_landscape(assembly$community)
    
    community_landscape <- community_landscape %>%
      dplyr::mutate(resident = ifelse(community_landscape$lma %in% community$lma, TRUE, FALSE)) 
    
    history_landscape <- purrr::map2(history_landscape, history, ~.x %>%
                                dplyr::mutate(resident = ifelse(.x$lma %in% .y$lma, TRUE, FALSE)))
  }
  
  
  tidied_assembly = list(history = history,
                         community = community,
                         history_landscape = history_landscape,
                         community_landscape = community_landscape)
  return(tidied_assembly)
}

tidy_community <- function(community){
  dplyr::tibble(community$traits %>% as_tibble(), 
         births = community$birth_rate,
         resident = TRUE)
}

tidy_history_community <- function(community_history){
  out <- purrr::map(community_history, tidy_community)
  
  out[[1]] %>%
    dplyr::select(-births, -resident) %>%
    names() -> traits

  
  for(i in seq_along(out)){
    if(i == 1){
      out[[i]] <- out[[i]] %>%
        dplyr::mutate(resident = FALSE)
    } else{
      trait_previous <- out[[i - 1]] %>%
        dplyr::unite(trait_whole, any_of(traits))
      
      trait_current <- out[[i]] %>%
        dplyr::unite(trait_whole, any_of(traits))
      
      out[[i]] <- out[[i]] %>%
        dplyr::mutate(resident = ifelse(trait_current$trait_whole %in% trait_previous$trait_whole,
                                TRUE, FALSE))
    }
    
  }
  
  
  out
}

