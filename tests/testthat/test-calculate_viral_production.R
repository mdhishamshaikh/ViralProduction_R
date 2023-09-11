# Config
expect_data_input <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('viralprod','data.frame'))
  expect_true(all(c('Location', 'Station_Number', 'Depth', 'Sample_Type', 'Timepoint', 'Replicate') %in% names(x)))
  expect_true(any(startsWith(names(x), 'c_')))
}

expect_output_vp_results <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared', 'VP_Method') %in% names(x)))
}

expect_output_vp_results_SR <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Time_Range', 'Population', 'Sample_Type', 'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared', 'VP_Method', 'Replicate') %in% names(x)))
}


# Perform
test_that("Wrapper function viral production calculation works", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', 
                                          package = "viralprod"))
  x <- new_viralprod_class(data_NJ2020_all)
  expect_data_input(x)
  
  calculate_viral_production(x = data_NJ2020_all, write_output = F) %>% expect_message()
  
  calculate_viral_production(x, write_output = F, SR_calc = F, BP_endpoint = F)
  expect_output_vp_results(.GlobalEnv$vp_results_output_df)
  
  data_NJ2020_less <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_less_populations.csv', 
                                           package = "viralprod"))
  x <- new_viralprod_class(data_NJ2020_less)
  expect_data_input(x)
  
  calculate_viral_production(x, write_output = F)
  expect_output_vp_results(.GlobalEnv$vp_results_output_df)
  expect_output_vp_results(.GlobalEnv$vp_results_output_T24_df)
  expect_output_vp_results_SR(.GlobalEnv$vp_results_output_SR_df)
  expect_output_vp_results(.GlobalEnv$vp_results_output_BP_df)
  
  calculate_viral_production(x, write_output = T, output_dir = '') %>% expect_error()
  
  data_NJ2020_more <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_more_populations.csv', 
                                           package = "viralprod"))
  x <- new_viralprod_class(data_NJ2020_more)
  expect_data_input(x)
  
  calculate_viral_production(x, write_output = F, SR_calc = F, BP_endpoint = F, methods = c(2,6,12))
  expect_output_vp_results(.GlobalEnv$vp_results_output_df)
  
  data_NJ2020_without_cViruses <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_without_cViruses.csv', 
                                                       package = "viralprod"))
  
  new_viralprod_class(data_NJ2020_without_cViruses) %>% expect_error()
})
