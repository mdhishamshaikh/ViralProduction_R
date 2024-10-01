#' Check populations to analyze
#' 
#' @description
#' Bacterial and viral counts are retrieved from flow cytometry data by selecting an area on the generated scatter plot,
#' a process called `gating`. During the gating process, different populations are defined based on the level of green
#' fluorescence and side scatter. Since gating is a manual process, the user is free to determine which populations to define. 
#' Depending on the area the gates encompass, count values are obtained and need to be as a column in the output data frame of the flow 
#' cytometry step named as followed: `c_PopulationName`. Given the output data frame of the flow cytometry step, 
#' the different populations to analyze, determined by the gating process, are defined. 
#' 
#' @param data Data frame with the output of the flow cytometry.
#'
#' @return A character vector with the different populations to analyze. An error occurs when the total virus population, `c_Viruses`, is not defined during the gating process. 
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
#' 
#' vp_check_populations(data_NJ2020_all)
#' 
#' # Case 2: Less populations defined during gating 
#' # (only the total viral and bacterial population for example)
#' data_NJ2020_less <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_less_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_less)
#' 
#' # Case 3: More populations defined during gating 
#' # (more niche populations, like c_V4, for example)
#' data_NJ2020_more <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_more_populations.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_more)
#' 
#' # Case 4: Total virus population is not defined during gating 
#' # (error expected)
#' data_NJ2020_without_cViruses <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_without_cViruses.csv', package = "viralprod"))
#' 
#' vp_check_populations(data_NJ2020_without_cViruses)
#' }
vp_check_populations <- function(data){
  if ('c_Viruses' %in% colnames(data)){
    .GlobalEnv$populations_to_analyze <- colnames(data)[grep("^c_", colnames(data))]
    print(paste("Following populations will be analyzed:", paste(.GlobalEnv$populations_to_analyze, collapse = ", ")))
  } else {
    stop('Total virus population, column c_Viruses, is not gated in output data frame of flow cytometry. Not able to perform viral production calculation!')
  }
}


#' Adding unique time ranges of the assay
#' 
#' Given a data frame that consists of column, `Timepoint`, representing the different sampling points of the assay,
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
#' 
#' vp_check_populations(data_NJ2020_all)
#' 
#' NJ2020_SR <- vp_separate_replicate_dataframe(data_NJ2020_all, add_timepoints = F)
#' 
#' vp_add_timepoints(NJ2020_SR)
#' 
#' NJ2020_AVG <- vp_average_replicate_dataframe(data_NJ2020_all, add_timepoints = F)
#' 
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
#' 
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

# New peak/valley function using pracma

vp_determine_peaks_with_se <- function(counts, sem) {
  # Load the necessary library
  library(pracma)
  
  # Step 1: Add a large constant to shift all values to positive range
  shift_constant <- abs(min(counts)) + 10e10  # Ensure all values are positive
  shifted_counts <- counts + shift_constant
  
  # Step 2: Identify peaks and valleys in the shifted data
  peak_indices <- findpeaks(shifted_counts)
  valley_indices <- findpeaks(-shifted_counts)  # Use negative to find valleys
  
  # Extract positions of peaks and valleys
  peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
  valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys
  
  # Initialize vectors for valid indices
  valid_peaks <- c()
  valid_valleys <- c()
  
  # Step 3: Pair each valley with the nearest peak to the right and validate with SEM
  for (valley_pos in valley_positions) {
    # Find the first peak that comes after the valley
    next_peak_index <- which(peak_positions > valley_pos)[1]
    
    # Skip if no peak is found to the right
    if (is.na(next_peak_index)) next
    
    # Get the position of this next peak
    peak_pos <- peak_positions[next_peak_index]
    
    # Calculate means and SEMs for the original (non-shifted) data
    valley_mean <- counts[valley_pos]  # Adjust position due to boundary
    valley_sem <- sem[valley_pos]
    
    peak_mean <- counts[peak_pos]  # Adjust position due to boundary
    peak_sem <- sem[peak_pos]
    
    # Step 4: Check SEM overlap condition
    if ((valley_mean + valley_sem) < (peak_mean - peak_sem)) {
      # If no overlap, store the valid indices
      valid_peaks <- c(valid_peaks, peak_pos - 1)
      valid_valleys <- c(valid_valleys, valley_pos - 1)
    }
  }
  
  # Return a list containing vectors of peak and valley indices
  return(list(peaks = valid_peaks, valleys = valid_valleys))
}

