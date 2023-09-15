#' Determine bacterial endpoint of assay
#' 
#' @description
#' In the VP samples, an increase of collision rates between the bacteriophages (viruses) and bacteria is noticed. 
#' This is probably due to the net increase in bacterial growth that is established during the assay. In VPC samples, 
#' this phenomenon isn't presented since treatment with antibiotic `Mitomycin-C` inhibits the
#' growth of bacteria. Such increased collision rates can result in an overestimated lytic viral production 
#' and a negative lysogenic viral production.
#' 
#' To mitigate this, we define the bacterial endpoint as the point at which the bacterial generation time becomes 
#' less than the assay duration, indication the appropriate moment to stop the assay and obtain less biased results.
#' The generation time is calculated based of the net increase of the total bacteria population in 
#' the VP samples of the assay.
#' 
#' @param data Data frame with the output of the flow cytometry.
#' @param visual If \code{FALSE}, a character with the time range to stop the assay is returned. The character value is from the form: T0_TX. If \code{TRUE}, an integer with the index of the time point to stop the assay is returned. (Default = \code{FALSE})
#' 
#' @return A character or integer value defining the time range to stop the assay to retrieve comparable results.
#' 
#' @name vp_bacterial_endpoint
#' @rdname vp_bacterial_endpoint
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_all)
#' 
#' # Bacterial endpoint is determined per station/experiment
#' subset_data <- subset(data_NJ2020_all, data_NJ2020_all$Station_Number == 2)
#' 
#' vp_bacterial_endpoint(subset_data)
#' vp_bacterial_endpoint(subset_data, visual = T)
#' }
vp_bacterial_endpoint <- function(data,
                                  visual = FALSE){
  DF_bacterial_and_VP_samples_only <- data %>%
    vp_average_replicate_dataframe(add_timepoints = F) %>%
    dplyr::filter(.data$Sample_Type == 'VP', .data$Microbe == 'Bacteria')
  
  bacterial_generation_time <- list()
  unique_timepoints <- unique(DF_bacterial_and_VP_samples_only$Timepoint)
  
  for (bacteria in unique(DF_bacterial_and_VP_samples_only$Population)){
    for (time in 2:length(unique_timepoints)){ # Calculate generation time relative towards starting point: time point 1 = T_0
      DF_current_bacteria <- DF_bacterial_and_VP_samples_only %>%
        dplyr::filter(.data$Population == bacteria) %>%
        dplyr::arrange(.data$Timepoint)
      
      generation_time <- (log10(2)*(DF_current_bacteria$Timepoint[time] - DF_current_bacteria$Timepoint[1])) / (log10(DF_current_bacteria$Mean[time]) - log10(DF_current_bacteria$Mean[1]))
      result <- c(unique_timepoints[time], bacteria, generation_time)
      bacterial_generation_time[[length(bacterial_generation_time) + 1]] <- result
    }
  }
  DF_bacterial_endpoint <- data.frame(t(sapply(bacterial_generation_time, c)))
  colnames(DF_bacterial_endpoint) <- c('Timepoint', 'Population', 'Generation_Time')
  
  DF_bacterial_endpoint <- DF_bacterial_endpoint %>%
    tidyr::pivot_wider(names_from = 'Population', values_from = 'Generation_Time')  %>%
    as.data.frame() %>%
    dplyr::mutate_at(dplyr::vars(-'Timepoint'), as.numeric)
  
  bacterial_endpoint_index <- intersect(which(DF_bacterial_endpoint$c_Bacteria > 0), which(DF_bacterial_endpoint$c_Bacteria < 24))[1]
  
  if (visual == T){
    if(is.na(bacterial_endpoint_index)){
      stop_assay <- length(unique_timepoints) - 1
    } else {
      stop_assay <- bacterial_endpoint_index - 1
    }
  } else {
    if (is.na(bacterial_endpoint_index)){
      stop_assay <- paste0("T", unique_timepoints[1], "_T", unique_timepoints[length(unique_timepoints)])
    } else {
      stop_assay <-  paste0("T", unique_timepoints[1], "_T", unique_timepoints[bacterial_endpoint_index])
    }
  }
  return(stop_assay)
}
