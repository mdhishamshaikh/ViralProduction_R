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
test_that("determine_vp_VIPCAL_separate_replicates", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
  expect_df_input_vp_VIPCAL(DF_SR)
  
  determine_vp_VIPCAL_separate_replicates(DF_SR) %>% expect_output_determine_vp_VIPCAL_separate_replicates()
  determine_vp_VIPCAL_separate_replicates(DF_SR) %>% ncol() %>% expect_equal(10)
})

test_that("determine_vp_VIPCAL_average_replicates", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_AVG <- vp_average_replicate_dataframe(data_NJ2020)
  expect_df_input_vp_VIPCAL(DF_AVG)
  
  determine_vp_VIPCAL_average_replicates(DF_AVG) %>% expect_output_determine_vp_VIPCAL_average_replicates()
  determine_vp_VIPCAL_average_replicates(DF_AVG) %>% ncol() %>% expect_equal(9)
})

test_that("determine_vp_VIPCAL_average_replicates_SE", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_AVG <- vp_average_replicate_dataframe(data_NJ2020)
  
  determine_vp_VIPCAL_average_replicates_SE(DF_AVG) %>% expect_output_determine_vp_VIPCAL_SE()
  determine_vp_VIPCAL_average_replicates_SE(DF_AVG) %>% ncol() %>% expect_equal(10)
})

test_that("determine_vp_VIPCAL_LMER_model", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
  
  determine_vp_VIPCAL_LMER_model(DF_SR) %>% expect_output_determine_vp_VIPCAL_average_replicates()
  determine_vp_VIPCAL_LMER_model(DF_SR) %>% ncol() %>% expect_equal(9)
})

test_that("determine_vp_VIPCAL_LMER_model_SE", {
  data_NJ2020 <- read.csv(system.file('extdata', 'NJ2020_subset.csv', package = "viralprod"))
  DF_SR <- vp_separate_replicate_dataframe(data_NJ2020)
  
  determine_vp_VIPCAL_LMER_model_SE(DF_SR) %>% expect_output_determine_vp_VIPCAL_SE()
  determine_vp_VIPCAL_LMER_model_SE(DF_SR) %>% ncol() %>% expect_equal(10)
})

