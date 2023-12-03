#' Tidies output from assembler_run.
#'
#' @param assembly 
#'
#' @return Returns the birth rate of the final community as well as the historical timeseries of birth rates for residents. If fitness landscapes are available, provides those as well.
#' @export
#'
#' @examples
tidy_assembly <- function(assembly){
  
  history <- tidy_history_community(assembly$history)
  community <- tidy_community(assembly$community)
  history_landscape <- NA
  community_landscape <- NA
  if(!is.null(assembly$community$fitness_points)){
    history_landscape <- tidy_landscape(assembly$history)
    community_landscape <- tidy_landscape(assembly$community)
    
    community_landscape <- community_landscape %>%
      mutate(resident = ifelse(community_landscape$lma %in% community$lma, TRUE, FALSE)) 
    
    history_landscape <- map2(history_landscape, history, ~.x %>%
                                mutate(resident = ifelse(.x$lma %in% .y$lma, TRUE, FALSE)))
  }
  
  
  tidied_assembly = list(history = history,
                         community = community,
                         history_landscape = history_landscape,
                         community_landscape = community_landscape)
  return(tidied_assembly)
}

tidy_community <- function(community){
  tibble(community$traits %>% as_tibble(), 
         births = community$birth_rate,
         invader = FALSE)
}

tidy_history_community <- function(community_history){
  out <- map(community_history, tidy_community)
  
  out[[1]] %>%
    select(-births, -invader) %>%
    names() -> traits

  
  for(i in seq_along(out)){
    if(i == 1){
      out[[i]] <- out[[i]] %>%
        mutate(invader = TRUE)
    } else{
      trait_previous <- out[[i - 1]] %>%
        unite(trait_whole, any_of(traits))
      
      trait_current <- out[[i]] %>%
        unite(trait_whole, any_of(traits))
      
      out[[i]] <- out[[i]] %>%
        mutate(invader = ifelse(trait_current$trait_whole %in% trait_previous$trait_whole,
                                FALSE, TRUE))
    }
    
  }
  
  
  out
}

tidy_landscape <- function(community){
  if(!is.null(names(community))){tidy_landscape_community(community) -> tidy_obj}
  if(is.null(names(community))){purrr::map(community, ~tidy_landscape_community(.x))-> tidy_obj}
  tidy_obj
}

tidy_landscape_community <- function(community){
  community$fitness_points %>%
    as_tibble()
}
