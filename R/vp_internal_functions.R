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
#' NJ2020 <- read.csv("./inst/NJ2020_subset.csv")
#' df <- df_SR(NJ2020, add_tp = F)
#' vp_add_timepoitns(df)
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

