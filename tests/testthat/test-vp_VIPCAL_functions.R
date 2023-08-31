# Config
expect_df_input_vp_VIPCAL <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true('Time_Range' %in% names(x))
}

expect_output_determine_vp_VIPCAL_separate_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Replicate', 'VP', 'abs_VP') %in% names(x)))
}

expect_output_determine_vp_VIPCAL_average_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('VP', 'abs_VP') %in% names(x)))
}

expect_output_determine_vp_VIPCAL_SE <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('VP', 'abs_VP', 'VP_SE') %in% names(x)))
}

# Perform
test_that("All VIPCAL functions work correctly", {
  data_NJ2020_all <- read.csv(system.file('extdata', 'NJ2020_Station_2_and_6_all_populations.csv', package = "viralprod"))
  vp_check_populations(data_NJ2020_all)
  
  DF_SR <- vp_separate_replicate_dataframe(data_NJ2020_all)
  expect_df_input_vp_VIPCAL(DF_SR)
  
  DF_AVG <- vp_average_replicate_dataframe(data_NJ2020_all)
  expect_df_input_vp_VIPCAL(DF_AVG)
  
  determine_vp_VIPCAL_separate_replicates(DF_SR) %>% expect_output_determine_vp_VIPCAL_separate_replicates()
  determine_vp_VIPCAL_separate_replicates(DF_SR) %>% ncol() %>% expect_equal(10)
  
  determine_vp_VIPCAL_average_replicates(DF_AVG) %>% expect_output_determine_vp_VIPCAL_average_replicates()
  determine_vp_VIPCAL_average_replicates(DF_AVG) %>% ncol() %>% expect_equal(9)
  
  determine_vp_VIPCAL_average_replicates_SE(DF_AVG) %>% expect_output_determine_vp_VIPCAL_SE()
  determine_vp_VIPCAL_average_replicates_SE(DF_AVG) %>% ncol() %>% expect_equal(10)
  
  determine_vp_VIPCAL_LMER_model(DF_SR) %>% expect_output_determine_vp_VIPCAL_average_replicates()
  determine_vp_VIPCAL_LMER_model(DF_SR) %>% ncol() %>% expect_equal(9)
  
  determine_vp_VIPCAL_LMER_model_SE(DF_SR) %>% expect_output_determine_vp_VIPCAL_SE()
  determine_vp_VIPCAL_LMER_model_SE(DF_SR) %>% ncol() %>% expect_equal(10)
})
