# Config
expect_df_input_vp_linear <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, c('data.table', 'data.frame'))
  expect_true('Time_Range' %in% names(x))
}

expect_output_determine_vp_linear_average_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('VP', 'abs_VP', 'VP_SE', 'VP_R_Squared') %in% names(x)))
}

expect_output_determine_vp_linear_separate_replicates <- function(x){
  expect_type(x, 'list')
  expect_s3_class(x, 'data.frame')
  expect_true(all(c('Replicate', 'VP', 'abs_VP', 'VP_SE', 'VP_R_Squared') %in% names(x)))
}

# Perform
test_that("Slope functions work that use separate replicate dataframe", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
  expect_df_input_vp_linear(DF_SR)
  
  determine_vp_linear_allpoints(DF_SR) %>% expect_output_determine_vp_linear_average_replicates()
  determine_vp_linear_allpoints(DF_SR) %>% ncol() %>% expect_equal(11)
  
  determine_vp_linear_separate_replicates(DF_SR) %>% expect_output_determine_vp_linear_separate_replicates()
  determine_vp_linear_separate_replicates(DF_SR) %>% ncol() %>% expect_equal(12)
  
  determine_vp_linear_LMER_model(DF_SR) %>% expect_output_determine_vp_linear_average_replicates()
  determine_vp_linear_LMER_model(DF_SR) %>% ncol() %>% expect_equal(11)
})

test_that("Slope functions work that use average replicate dataframe", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_AVG <- vp_average_replicate_dataframe(data_NJ2020)
  expect_df_input_vp_linear(DF_AVG)
  
  determine_vp_linear_average_replicates(DF_AVG) %>% expect_output_determine_vp_linear_average_replicates()
  determine_vp_linear_average_replicates(DF_AVG) %>% ncol() %>% expect_equal(11)
})
