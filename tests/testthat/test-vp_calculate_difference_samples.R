# Config
expect_df_input_vp_calc_diff_linear <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('tag', 'Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type',
                    'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared') %in% names(x)))
}

expect_df_input_vp_calc_diff_VIPCAL <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('tag', 'Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type',
                    'VP', 'abs_VP') %in% names(x)))
}

expect_df_input_vp_calc_diff_VIPCAL_SE <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('tag', 'Location', 'Station_Number', 'Depth', 'Time_Range', 'Population', 'Sample_Type',
                    'VP', 'abs_VP', 'VP_SE') %in% names(x)))
}

# Perform
test_that("Difference samples can be calculated for all methods", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
  vp_check_populations(data_NJ2020_all)
  
  DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
  DF_AVG <- vp_average_replicate_dataframe(data_NJ2020_all) %>% subset(Sample_Type != 'Diff')
  
  vp_linear <- determine_vp_linear_allpoints(DF_SR)
  expect_df_input_vp_calc_diff_linear(vp_linear)
  
  test_linear <- vp_calculate_difference_samples(vp_linear)
  expect_true(any(test_linear$Sample_Type == "Diff"))
  
  vp_VIPCAL <- determine_vp_VIPCAL_average_replicates(DF_AVG)
  expect_df_input_vp_calc_diff_VIPCAL(vp_VIPCAL)
  
  test_VIPCAL <- vp_calculate_difference_samples(vp_VIPCAL, VIPCAL = T)
  expect_true(any(test_VIPCAL$Sample_Type == "Diff"))
  
  vp_VIPCAL_SE <- determine_vp_VIPCAL_average_replicates_SE(DF_AVG)
  expect_df_input_vp_calc_diff_VIPCAL_SE(vp_VIPCAL_SE)
  
  test_VIPCAL_SE <- vp_calculate_difference_samples(vp_VIPCAL_SE, VIPCAL = T, SE = T)
  expect_true(any(test_VIPCAL_SE$Sample_Type == "Diff"))
})
