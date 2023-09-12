#' Check populations to analyze
#' 
#' @description
#' Bacterial and viral counts are retrieved from flow cytometry data by selecting an area on the generated scatter plot,
#' this process is called `gating`. During the gating process, different populations are defined based on the side scatter
#' and green fluorescence. Since gating is a manual process, the user is free to determine which populations to define. 
#' The count values, retrieved from the scatter plot, need to be as a column in the output data frame of the flow 
#' cytometry step named as followed: `c_PopulationName`. Given the output data frame of the flow cytometry step, 
#' the different populations to analyze, determined by the gating process, are defined. 
#' 
#' @param data Data frame with the output of the flow cytometry.
#'
#' @return A character vector with the different populations to analyze will be available in the global environment. An error occurs when the total virus population, `c_Viruses`, is not defined during the gating process. 
#' @export
#' 
#' @name vp_check_populations
#' @rdname vp_check_populations
#'
#' @examples \dontrun{
#' # Case 1: General, 7 different populations defined during gating 
#' # (most common populations)
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_all)
#' 
#' # Case 2: Less populations defined during gating 
#' # (only the total viral and bacterial population for example)
#' data_NJ2020_less <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_less_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_less)
#' 
#' # Case 3: More populations defined during gating 
#' # (more niche populations, like c_V4, for example)
#' data_NJ2020_more <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_more_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_more)
#' 
#' # Case 4: Total virus population is not defined during gating 
#' # (error expected)
#' data_NJ2020_without_cViruses <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_without_cViruses.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_without_cViruses)
#' }
vp_check_populations <- function(data){
  if ('c_Viruses' %in% colnames(data)){
    .GlobalEnv$populations_to_analyze <- colnames(data)[grep("^c_", colnames(data))]
    message(paste("Following populations will be analyzed:", paste(.GlobalEnv$populations_to_analyze, collapse = ", ")))
  } else {
    stop('Total virus population, column c_Viruses, is not gated in output data frame of flow cytometry. Not able to perform viral production calculation!')
  }
}


#' Adding unique time ranges of the assay
#' 
#' Given a data frame that consists of column, `Timepoint`, that represents the different sampling points of the assay,
#' a column with the different time ranges of the assay is added to the original data frame. 
#'
#' @param DF Data frame with the count for each population and each sample at the different time points of the assay.
#'
#' @return Expanded data frame with time ranges added as new column.
#' 
#' @name vp_add_timepoints
#' @rdname vp_add_timepoints
#' @noRd
#'
#' @examples 
#' x <- data.frame(Timepoint = c(0,3,6,9,12,24))
#' vp_add_timepoints(x)
#' 
#' \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_all)
#' 
#' NJ2020_SR <- vp_separate_replicate_dataframe(data_NJ2020_all, add_timepoints = F)
#' vp_add_timepoints(NJ2020_SR)
#' 
#' NJ2020_AVG <- vp_average_replicate_dataframe(data_NJ2020_all, add_timepoints = F)
#' vp_add_timepoints(NJ2020_AVG)
#' }
vp_add_timepoints <- function(DF){
  timepoints <- unique(as.numeric(DF$Timepoint))
  
  colnames<- c() 
  for(col in 2:length(timepoints)){
    timerange_name <- paste("T", timepoints[1], "_T", timepoints[col], sep = "")
    colnames[length(colnames)+1] <- timerange_name
  }
  
  colvalues<- c() 
  for(col in 2:length(timepoints)){
    timerange_value <- paste("T", timepoints[1], ":T", timepoints[col], sep = "")
    colvalues[length(colvalues)+1] <- timerange_value
  }
  
  number_of_columns <- ncol(DF)
  DF[colnames] <- NA
  
  for (i in 1:(length(timepoints)-1)) { 
    conditions <- DF$Timepoint %in% timepoints[1:(i+1)]
    DF[conditions, (number_of_columns + i)] <- colvalues[i]
  }
  
  DF <- DF %>%
    tidyr::pivot_longer(cols = dplyr::all_of(colnames), names_to = "Time_Range", values_to = "Time_Time") %>%
    tidyr::drop_na()
  
  return(DF)
}


#' Determine peaks and valleys
#' 
#' @description
#' `VIPCAL` calculates the viral production based on the average of increments. To get these increments, peaks
#' and valleys in the viral counts need to be determined. VIPCAL has his own issues, namely that the standard error 
#' has a big influence on the results. `VIPCAL-SE` goes one step further and takes the standard error into account
#' when determining peaks and valleys. Because of that, only TRUE increments (increments without overlapping standard errors) are returned. 
#'
#' @param count_values Column of data frame with viral count values.
#' @param count_se Column of data frame with standard error on the viral count values.
#'
#' @return Vector with the indices of the peaks or valleys in the count data.
#' 
#' @name vp_peaks_and_valleys
#' @rdname vp_peaks_and_valleys
#' @noRd
#'
#' @examples \dontrun{
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' vp_check_populations(data_NJ2020_all)
#' 
#' DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
#' DF_AVG <- vp_average_replicate_dataframe(data_NJ2020_all)
#' 
#' # Adding two values to make sure the first and last element of the count values are not dismissed
#' vp_determine_peaks(c(+10e+10, DF_SR$Count, -10e+10))
#' vp_determine_peaks(c(+10e+10, DF_AVG$Mean, -10e+10))
#' 
#' vp_determine_valleys(c(+10e+10, DF_SR$Count, -10e+10))
#' vp_determine_valleys(c(+10e+10, DF_AVG$Mean, -10e+10))
#' 
#' vp_determine_peaks_with_se(c(+10e+10, DF_AVG$Mean, -10e+10),
#'                            c(0, DF_AVG$SE, 0))
#' vp_determine_valleys_with_se(c(+10e+10, DF_AVG$Mean, -10e+10),
#'                              c(0, DF_AVG$SE, 0))
#' }
vp_determine_peaks <- function(count_values){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){ 
    sign_index <- sign(count_values[index+1] - count_values[index])
    result_list[length(result_list)+1] <- sign_index 
  }
  
  return(which(diff(result_list) < 0)) # PEAK if difference is negative
}


#' @rdname vp_peaks_and_valleys
#' @noRd
vp_determine_valleys <- function(count_values){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){ 
    sign_index <- sign(count_values[index+1] - count_values[index]) 
    result_list[length(result_list)+1] <- sign_index 
  }
  
  return(which(diff(result_list) > 0)) # VALLEY if difference is positive
}

#' @rdname vp_peaks_and_valleys
#' @noRd
vp_determine_peaks_with_se <- function(count_values, 
                                       count_se){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){
    sign_index <- sign((count_values[index+1] - count_se[index+1]) - (count_values[index] + count_se[index]))
    result_list[length(result_list)+1] <- sign_index
  }
  return(which(diff(result_list) < 0)) # PEAK if difference is negative
}


#' @rdname vp_peaks_and_valleys
#' @noRd
vp_determine_valleys_with_se <- function(count_values, 
                                         count_se){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){
    sign_index <- sign((count_values[index+1] - count_se[index+1]) - (count_values[index] + count_se[index]))
    result_list[length(result_list)+1] <- sign_index
  }
  return(which(diff(result_list) > 0)) # VALLEY if difference is positive
}