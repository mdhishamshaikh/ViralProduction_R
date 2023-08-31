#' Construct viral count data frame
#' 
#' @description
#' Given the output data frame of the flow cytometry step, produce the correct data frame for further
#' calculations. Based on the replicate treatment, a different data frame will be produced. A separate
#' replicate treatment will take each of the replicates into account, as a average replicate treatment
#' will average over the replicates to retrieve an average count value. 
#' 
#' `vp_separate_replicate_dataframe` creates data frame by taking the separate replicates into account.
#' 
#' `vp_average_replicate_dataframe` creates data frame by taking the average over the replicates.
#' 
#' @param data Data frame with the output of the flow cytometry.
#' @param keep_0.22_samples If \code{FALSE}, 0.22 samples will be removed from data. These represent the control samples and are normally always filtered out. If you want to keep control samples, set to \code{TRUE}. (Default = \code{FALSE})
#' @param add_timepoints If \code{TRUE}, different time ranges will be added with \code{vp_add_timepoints()}. If \code{FALSE}, no addition of time ranges to data frame. (Default = \code{TRUE})
#'
#' @return Data frame with the count for each population and each sample at the different time points of the assay. If \code{add_timepoints = T}, a column with the time ranges is added to the data frame. 
#' 
#' @name vp_dataframes
#' @rdname vp_dataframes
#' 
#' @importFrom rlang !!!
#' 
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_all)
#' 
#' vp_separate_replicate_dataframe(data_NJ2020_all)
#' vp_separate_replicate_dataframe(data_NJ2020_all, add_timepoints = F)
#' vp_separate_replicate_dataframe(data_NJ2020_all, keep_0.22_samples = T)
#' 
#' vp_average_replicate_dataframe(data_NJ2020_all)
#' vp_average_replicate_dataframe(data_NJ2020_all, add_timepoints = F)
#' }
vp_separate_replicate_dataframe <- function(data, keep_0.22_samples = FALSE, add_timepoints = TRUE){
  SR_dataframe <- data %>%
    dplyr::select(dplyr::all_of(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', !!!.GlobalEnv$populations_to_analyze))) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with('c_'), names_to = 'Population', values_to = 'Count') %>% 
    dplyr::mutate(Microbe = dplyr::if_else(.data$Population %in% c('c_Bacteria', 'c_HNA', 'c_LNA'), 'Bacteria', 'Viruses')) %>%
    tidyr::unite(dplyr::all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F) %>%
    dplyr::arrange('tag', 'Sample_Type','Replicate','Population', as.numeric(.data$Timepoint))
  
  if (keep_0.22_samples == F){
  SR_dataframe <- SR_dataframe[SR_dataframe$Sample_Type != '0.22',]
  }
  
  if (add_timepoints == T){
    result_list <- list()
    
    for(combi_tag in unique(SR_dataframe$tag)){
      SR_dataframe_2 <- SR_dataframe %>%
        dplyr::filter(.data$tag == combi_tag)
      
      SR_dataframe_2 <- vp_add_timepoints(SR_dataframe_2)
      result_list[[length(result_list)+1]] <- SR_dataframe_2
    }
    
    SR_dataframe_with_timepoints <- data.table::rbindlist(result_list)
  }else {
    SR_dataframe_with_timepoints <- as.data.frame(SR_dataframe)
  }
  return(SR_dataframe_with_timepoints)
}


#' @rdname vp_dataframes
vp_average_replicate_dataframe <- function(data, add_timepoints = TRUE){
  dataframe_without_controls <- data[data$Sample_Type != '0.22',]
  
  AVG_dataframe <- dataframe_without_controls %>%
    dplyr::select(dplyr::all_of(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', !!!.GlobalEnv$populations_to_analyze))) %>%
    tidyr::pivot_longer(cols = dplyr::starts_with('c_'), names_to = 'Population', values_to = 'Count') %>%
    dplyr::group_by(.data$Location, .data$Station_Number, .data$Depth, .data$Sample_Type, .data$Timepoint, .data$Population) %>%
    dplyr::summarise(n = dplyr::n(), Mean = mean(.data$Count), SE = plotrix::std.error(.data$Count))
  
  AVG_dataframe_only_means <- AVG_dataframe %>%
    dplyr::select(-'SE') %>%
    tidyr::spread('Sample_Type', 'Mean')
  
  AVG_dataframe_only_se <- AVG_dataframe %>%
    dplyr::select(-'Mean') %>%
    tidyr::spread('Sample_Type', 'SE')
  
  if ('VPC' %in% AVG_dataframe$Sample_Type){
    AVG_dataframe_only_means$Diff <- with(AVG_dataframe_only_means, VPC - VP)
    AVG_dataframe_only_means <- tidyr::pivot_longer(AVG_dataframe_only_means, cols = dplyr::all_of(c('VP', 'VPC', 'Diff')), names_to = 'Sample_Type', values_to = 'Mean')
    
    AVG_dataframe_only_se$Diff <- with(AVG_dataframe_only_se, VPC + VP)
    AVG_dataframe_only_se <- tidyr::pivot_longer(AVG_dataframe_only_se, cols = dplyr::all_of(c('VP', 'VPC', 'Diff')), names_to = 'Sample_Type', values_to = 'SE')
  }
  
  AVG_dataframe_merged <- merge(AVG_dataframe_only_means, AVG_dataframe_only_se, 
                                by = c('Location', 'Station_Number', 'Depth', 'Timepoint', 'Population', 'n', 'Sample_Type')) %>%
    dplyr::mutate(Microbe = dplyr::if_else(.data$Population %in% c('c_Bacteria', 'c_HNA', 'c_LNA'), 'Bacteria', 'Viruses')) %>%
    tidyr::unite(dplyr::all_of(c('Location', 'Station_Number', 'Depth')), col = 'tag', remove = F) %>%
    dplyr::arrange('tag', 'Sample_Type','Replicate','Population', as.numeric(.data$Timepoint))
  
  if (add_timepoints == TRUE){
    result_list <- list()
    
    for(combi_tag in unique(AVG_dataframe_merged$tag)){
      AVG_dataframe_merged_2 <- AVG_dataframe_merged %>%
        dplyr::filter(.data$tag == combi_tag)
      
      AVG_dataframe_merged_2 <- vp_add_timepoints(AVG_dataframe_merged_2)
      result_list [[length(result_list)+1]] <- AVG_dataframe_merged_2
    }
    AVG_dataframe_with_timepoints <- data.table::rbindlist(result_list)
  }else { 
    AVG_dataframe_with_timepoints <- as.data.frame(AVG_dataframe_merged)
  }
  
  return(AVG_dataframe_with_timepoints)
}
