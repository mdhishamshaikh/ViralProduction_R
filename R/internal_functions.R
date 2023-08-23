#' Adding time ranges of assay to dataframe
#'
#' @param DF Dataframe consisting of a numerical column Timepoint that has the different timepoints
#'
#' @return The input dataframe added with two columns Time_Range and Time_Time that show the different time ranges of the assay
#' @export
#'
#' @examples
#' x <- data.frame(Timepoint = c(0,3,6,9,12,24))
#' tp(x)
tp <- function(DF){
  
  TP <- unique(as.numeric(DF$Timepoint)) # Unique timepoints of assay
  
  # Timepoint always starts from T = 0 => iteration starts from 2nd timepoint
  colnames<- c() # Of the form: "TX_TY"
  
  for(col in 2:length(TP)){
    a <- paste("T", TP[1], "_T", TP[col], sep = "") # Make current timepoint
    colnames[length(colnames)+1] <- a # Add it as column to vector colnames
  }
  
  colvalues<- c() # Of the form: "TX:TY"
  for(col in 2:length(TP)){
    a <- paste("T", TP[1], ":T", TP[col], sep = "") # Make current timepoint
    colvalues[length(colvalues)+1] <- a # Add it as column to vector colvalues
  }
  
  # Adding columns to link colnames to colvalues
  ncol <- ncol(DF)
  DF[colnames] <- NA # Initiate new columns
  
  # Simplified loop for adding the columns to dataframe
  for (i in 1:(length(TP)-1)) { 
    conditions <- DF$Timepoint %in% TP[1:(i+1)]
    DF[conditions, (ncol + i)] <- colvalues[i]
  }
  
  # We have added columns to each replicate with its timeframe
  # Change the format => increase number of rows and decrease the number of columns with pivot_longer()
  DF <- DF %>%
    tidyr::pivot_longer(cols = dplyr::all_of(colnames), names_to = "Time_Range", values_to = "Time_Time") %>%
    tidyr::drop_na() # Drops all the rows which have a NA value
  
  return(DF)
}