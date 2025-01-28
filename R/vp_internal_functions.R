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
#' vp_determine_peaks_valleys_pracma(c(+10e+100, DF_SR$Count+10e+10, -10e+100))[[1]]
#' vp_determine_peaks_valleys_pracma(c(+10e+100, DF_AVG$Mean+10e+10, -10e+100))[[2]]
#' 
#' vp_determine_peaks_valleys_pracma(c(+10e+100, DF_SR$Count+10e+10, -10e+100))[[1]]
#' vp_determine_peaks_valleys_pracma(c(+10e+100, DF_AVG$Mean+10e+10, -10e+100))[[2]]
#' 
#' vp_determine_peaks_valleys_pracma(c(+10e+100, DF_AVG$Mean+10e+10, -10e+100),
#'                            c(1, DF_AVG$SE, 1))[[1]]
#' vp_determine_peaks_valleys_pracma(c(+10e+100, DF_AVG$Mean+10e+10, -10e+100),
#'                              c(1, DF_AVG$SE, 1))[[2]]
#' }

# New peak/valley function using pracma

vp_determine_peaks_valleys_pracma <- function(counts) {
 
  # Step 1: Identifying initial peaks and valleys in the data
  peak_indices <- pracma::findpeaks(counts)
  valley_indices <- pracma::findpeaks(-counts)  # Using negative to find valleys
  
  # Extracting positions of peaks and valleys
  peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
  valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys
  
  # Initializing vectors for the paired indices
  paired_peaks <- c()
  paired_valleys <- c()
  
  # Step 2: Pairing each valley with the nearest peak to the right
  for (valley_pos in valley_positions) {
    # Findig the first peak that comes after the valley
    next_peak_index <- which(peak_positions > valley_pos)[1]
    
    # Skipping if no peak is found to the right
    if (is.na(next_peak_index)) next
    
    # Getting the position of this next peak and add to the list
    peak_pos <- peak_positions[next_peak_index]
    paired_peaks <- c(paired_peaks, peak_pos)
    paired_valleys <- c(paired_valleys, valley_pos)
  }
  
  # Correcting for flanking constants
  paired_peaks <- paired_peaks -1
  paired_valleys <- paired_valleys -1
  
  # Returning a list containing vectors of paired peak and valley indices
  return(list(peaks = paired_peaks, valleys = paired_valleys))
}


#' @rdname vp_peaks_and_valleys
#' @noRd


vp_determine_peaks_and_valleys_with_se_pracma <- function(counts, sem) {
  # Loading the necessary library
  library(pracma)
  # Step 0: Overwriting NAs in SEM to 1
  sem <- ifelse(is.na(sem), 1, sem)
  
  # Step 1: Identifying initial peaks and valleys in the data
  peak_indices <- pracma::findpeaks(counts)
  valley_indices <- pracma::findpeaks(-counts)  # Using negative to find valleys
  
  # Extractg positions of peaks and valleys
  peak_positions <- peak_indices[, 2]  # 2nd column contains the positions of peaks
  valley_positions <- valley_indices[, 2]  # 2nd column contains the positions of valleys
  
  # Initializing vectors for valid indices
  valid_peaks <- c()
  valid_valleys <- c()
  
  # Step 2: Pairng each valley with the nearest peak to the right and validate with SEM
  for (valley_pos in valley_positions) {
    # Find the first peak that comes after the valley
    next_peak_index <- which(peak_positions > valley_pos)[1]
    
    # Skipping if no peak is found to the right
    if (is.na(next_peak_index)) next
    
    # Getting the position of this next peak
    peak_pos <- peak_positions[next_peak_index]
    
    # Calculating means and SEMs for the original data
    valley_mean <- counts[valley_pos]
    valley_sem <- sem[valley_pos]
    peak_mean <- counts[peak_pos]
    peak_sem <- sem[peak_pos]
    
    # Checking SEM overlap condition
    if ((valley_mean + valley_sem) < (peak_mean - peak_sem)) {
      # If no overlap, store initial valid valley and peak indices
      valid_peaks <- c(valid_peaks, peak_pos)
      valid_valleys <- c(valid_valleys, valley_pos)
    }
  }
  
  # Step 3: Refinig valid peaks by moving left to find the last insignificant neighbor,
  # only up to the position just before the associated valley
  final_peaks <- c()
  for (i in seq_along(valid_peaks)) {
    peak_pos <- valid_peaks[i]
    valley_pos <- valid_valleys[i]
    current_peak <- peak_pos
    left_neighbor <- current_peak - 1
    
    # Moving leftward until a significant difference is found or we reach just before the valley
    while (left_neighbor > 0 && left_neighbor > valley_pos) {  # Ensure we don't move past the valley
      left_mean <- counts[left_neighbor]
      left_sem <- sem[left_neighbor]
      peak_mean <- counts[current_peak]
      peak_sem <- sem[current_peak]
      
      # Checking if current peak is insignificantly different from left neighbor
      if ((left_mean + left_sem) < (peak_mean - peak_sem)) {
        # Stopping if significantly different
        break
      }
      
      # Moving to left neighbor if insignificantly different
      current_peak <- left_neighbor
      left_neighbor <- current_peak - 1
    }
    
    # Storing the final validated peak position
    final_peaks <- c(final_peaks, current_peak)
  }
  
  # Step 4: Refinng valid valleys by moving right to find the last insignificant neighbor,
  # only up to the position just before the associated (refined) peak
  final_valleys <- c()
  for (i in seq_along(valid_valleys)) {
    valley_pos <- valid_valleys[i]
    current_peak <- final_peaks[i]  # Use the dynamically refined peak position for this valley
    current_valley <- valley_pos
    right_neighbor <- current_valley + 1
    
    # Moving rightward until a significant difference is found or we reach the refined peak position
    while (right_neighbor <= length(counts) && right_neighbor < current_peak) {  # Ensurinng valley does not reach the refined peak
      right_mean <- counts[right_neighbor]
      right_sem <- sem[right_neighbor]
      valley_mean <- counts[current_valley]
      valley_sem <- sem[current_valley]
      
      # Checking if current valley is insignificantly different from right neighbor
      if ((valley_mean + valley_sem) < (right_mean - right_sem)) {
        # Stop if significantly different
        break
      }
      
      # Moving to right neighbor if insignificantly different
      current_valley <- right_neighbor
      right_neighbor <- current_valley + 1
    }
    
    # Storing the final validated valley position
    final_valleys <- c(final_valleys, current_valley) 
  }
  
  # Adjusting for boundary additions if necessary
  final_peaks <- final_peaks - 1
  final_valleys <- final_valleys - 1
  
  # Returning a list containing vectors of refined peak and valley indices
  return(list(peaks = final_peaks, valleys = final_valleys))
}
#' @rdname vp_peaks_and_valleys
#' @noRd



