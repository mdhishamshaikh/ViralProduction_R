#' Adding timepoints
#' 
#' Given a dataframe that consists of column, 'Timepoint', that represents the different points of measuring the assay.
#' \code{\link{vp_add_timepoints()}} will add the different time ranges to the dataframe. 
#'
#' @param DF Dataframe.
#'
#' @return Same dataframe with time ranges added as new column.
#' 
#' @noRd
#'
#' @examples 
#' x <- data.frame(Timepoint = c(0,3,6,9,12,24))
#' vp_add_timepoints(x)
#' 
#' \dontrun{
#' NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' NJ2020_SR <- vp_separate_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_add_timepoints(NJ2020_SR)
#' 
#' NJ2020_AVG <- vp_average_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_add_timepoints(NJ2020_AVG)
#' }
vp_add_timepoints <- function(DF){
  timepoints <- unique(as.numeric(DF$Timepoint))
  
  # Determine time ranges
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
  
  # Adding to dataframe
  ncol <- ncol(DF)
  DF[colnames] <- NA
  
  for (i in 1:(length(timepoints)-1)) { 
    conditions <- DF$Timepoint %in% timepoints[1:(i+1)]
    DF[conditions, (ncol + i)] <- colvalues[i]
  }
  
  DF <- DF %>%
    tidyr::pivot_longer(cols = dplyr::all_of(colnames), names_to = "Time_Range", values_to = "Time_Time") %>%
    tidyr::drop_na() # Drops all the rows which have a NA value
  
  return(DF)
}


#' Determine peaks
#' 
#' VIPCAL calculates the viral production based off the average of increments. To get these increments peaks
#' and valleys in the viral counts need to be determined. \code{\link{vp_determine_peaks()}} takes count values
#' as input and determines the indexes of peaks.  
#'
#' @param count_values Column with viral count values.
#'
#' @return Vector with the indices of the peaks in the count data.
#' 
#' @noRd
#'
#' @examples \dontrun{
#' NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' # Adding two values to make sure the first and last element of the count values are not dismissed
#' NJ2020_SR <- vp_separate_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_determine_peaks(c(+10e+10, NJ2020_SR$Count, -10e+10))
#' 
#' NJ2020_AVG <- vp_average_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_determine_peaks(c(+10e+10, NJ2020_AVG$Mean, -10e+10))
#' }
vp_determine_peaks <- function(count_values){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){ 
    sign_index <- sign(count_values[index+1] - count_values[index])
    result_list[length(result_list)+1] <- sign_index 
  }
  
  return(which(diff(result_list) < 0)) # PEAK if difference is negative
}


#' Determine valleys
#' 
#' VIPCAL calculates the viral production based off the average of increments. To get these increments peaks
#' and valleys in the viral counts need to be determined. \code{\link{vp_determine_valleys()}} takes count values
#' as input and determines the indexes of valleys.  
#'
#' @param count_values Column with viral count values.
#'
#' @return Vector with the indices of the valleys in the count data.
#' 
#' @noRd
#'
#' @examples \dontrun{
#' NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' # Adding two values to make sure the first and last element of the count values are not dismissed
#' NJ2020_SR <- vp_separate_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_determine_valleys(c(+10e+10, NJ2020_SR$Count, -10e+10))
#' 
#' NJ2020_AVG <- vp_average_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_determine_valleys(c(+10e+10, NJ2020_AVG$Mean, -10e+10))
#' }
vp_determine_valleys <- function(count_values){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){ 
    sign_index <- sign(count_values[index+1] - count_values[index]) 
    result_list[length(result_list)+1] <- sign_index 
  }
  
  return(which(diff(result_list) > 0)) # VALLEY if difference is positive
}

#' Determine peaks with standard error taken into account
#' 
#' VIPCAL calculates the viral production based off the average of increments. To get these increments peaks
#' and valleys in the viral counts need to be determined. VIPCAL has his own issues, the standard error has a big
#' influence on the results. Some of the VIPCAL variants take the standard error into account when determining peaks
#' and valleys. Because of that, only TRUE increments (increments without overlapping standard errors) are returned. 
#' \code{\link{vp_determine_peaks_with_se()}} takes count values and corresponding standard errors as input and determines 
#' the indexes of peaks.  
#'
#' @param count_values Column with viral count values.
#' @param count_se Column with corresponding standard error of viral count value.
#'
#' @return Vector with the indices of the peaks in the count data.
#' 
#' @noRd
#'
#' @examples \dontrun{
#' NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' # Adding two values to make sure the first and last element of the count values are not dismissed
#' NJ2020_AVG <- vp_average_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_determine_peaks_with_se(c(+10e+10, NJ2020_AVG$Mean, -10e+10),
#'                            c(0, NJ2020_AVG$SE, 0))
#' }
vp_determine_peaks_with_se <- function(count_values, count_se){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){
    sign_index <- sign((count_values[index+1] - count_se[index+1]) - (count_values[index] + count_se[index]))
    result_list[length(result_list)+1] <- sign_index
  }
  return(which(diff(result_list) < 0)) # PEAK if difference is negative
}


#' Determine valleys with standard error taken into account
#' 
#' VIPCAL calculates the viral production based off the average of increments. To get these increments peaks
#' and valleys in the viral counts need to be determined. VIPCAL has his own issues, the standard error has a big
#' influence on the results. Some of the VIPCAL variants take the standard error into account when determining peaks
#' and valleys. Because of that, only TRUE increments (increments without overlapping standard errors) are returned. 
#' \code{\link{vp_determine_valleys_with_se()}} takes count values and corresponding standard errors as input and determines 
#' the indexes of valleys.  
#'
#' @param count_values Column with viral count values.
#' @param count_se Column with corresponding standard error of viral count value.
#'
#' @return Vector with the indices of the valleys in the count data.
#' 
#' @noRd
#'
#' @examples \dontrun{
#' NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
#' 
#' # Adding two values to make sure the first and last element of the count values are not dismissed
#' NJ2020_AVG <- vp_average_replicate_dataframe(NJ2020, add_timepoints = F)
#' vp_determine_valleys_with_se(c(+10e+10, NJ2020_AVG$Mean, -10e+10),
#'                            c(0, NJ2020_AVG$SE, 0))
#' }
vp_determine_valleys_with_se <- function(count_values, count_se){
  result_list <- c()
  
  for (index in 1:(length(count_values)-1)){
    sign_index <- sign((count_values[index+1] - count_se[index+1]) - (count_values[index] + count_se[index]))
    result_list[length(result_list)+1] <- sign_index
  }
  return(which(diff(result_list) > 0)) # VALLEY if difference is positive
}