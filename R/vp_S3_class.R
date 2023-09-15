#' Create S3 class for input data frames
#' 
#' @description
#' The `viralprod` package is dependent on two data frames, the output data frame of the flow cytometry step and a data 
#' frame with the original abundances in the sample. To make the R package more organized, maintainable and user-friendly, 
#' `S3 object-oriented programming` is introduced. Based on the requirements of the data frames, a separate S3 class, unique to this package, is defined.
#' The functions to calculate, analyze and visualize the viral production will only work on these specific
#' S3 classes, otherwise an error will be thrown. 
#' 
#' `vp_class_count_data` checks the output data frame of the flow cytometry step. If the requirements are met, `viralprod` class is added.
#' 
#' `vp_class_ori_abu` checks the original abundances data frame. If the requirements are met, `viralprod_analyze` class is added. 
#' 
#' @param data Data frame with the output of the flow cytometry.
#' @param original_abundances Data frame with the abundances of bacterial and virus population in the original sample.
#'
#' @return If data frame meets the requirements, a new class is added. Otherwise, an error will be thrown.
#' @export
#' 
#' @name vp_S3_class
#' @rdname vp_S3_class
#'
#' @examples \dontrun{
#' # Setup
#' data_NJ2020_all <- read.csv(system.file('extdata', 
#' 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
#' 
#' original_abundances_NJ2020 <- read.csv(system.file('extdata',
#' 'NJ2020_original_abundances.csv', package = "viralprod"))
#' 
#' # Perform
#' data_S3_class <- vp_class_count_data(data_NJ2020_all)
#' 
#' original_abundances_S3_class <- vp_class_ori_abu(original_abunances_NJ2020)
#' }
vp_class_count_data <- function(data = data.frame()){
  ## Assertions on the class
  # 1. Input data needs to be data frame
  checkmate::assert_data_frame(data)
  
  # 2. Data frame has to have certain columns
  vp_check_populations(data)
  has_correct_columns <- all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate', .GlobalEnv$populations_to_analyze) %in% names(data))
  checkmate::assert_true(has_correct_columns)
  
  # 3. Columns need to be of correct class
  if (has_correct_columns){
    checkmate::assert_true(all(sapply(data[c('Timepoint', .GlobalEnv$populations_to_analyze)], is.numeric)) && 
                             all(sapply(data[c('Station_Number', 'Depth', 'Replicate')], is.integer)) && 
                             all(sapply(data[c('Location', 'Sample_Type')], is.character)))
  }
  
  ## Return object with new class
  attr(data, "class") <- unique(c("viralprod", class(data)))
  
  return(data)
}


#' @export
#' @rdname vp_S3_class
vp_class_ori_abu <- function(original_abundances = data.frame()){
  ## Assertions on the class
  # 1. Input data needs to be data frame
  checkmate::assert_data_frame(original_abundances)
  
  # 2. Data frame has to have certain columns
  has_correct_columns <- all(c('Station_Number', 'Total_Bacteria', 'Total_Viruses') %in% names(original_abundances))
  checkmate::assert_true(has_correct_columns)
  
  # 3. Columns need to be of correct class
  if (has_correct_columns){
    checkmate::assert_true(all(sapply(original_abundances[c('Total_Bacteria', 'Total_Viruses')], is.numeric)) && 
                             all(sapply(original_abundances[c('Station_Number')], is.integer)))
  }
  
  ## Return object with new class
  attr(original_abundances, "class") <- unique(c("viralprod_analyze", class(original_abundances)))
  
  return(original_abundances)
}
